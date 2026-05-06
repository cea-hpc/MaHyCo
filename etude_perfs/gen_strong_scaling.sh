#!/usr/bin/env bash
# ==============================================================================
# gen_strong_scaling.sh
#
# Genere des fichiers .arc pour une etude de scalabilite forte.
# La balise <nsd>nx ny nz</nsd> est mise a jour pour chaque point de mesure.
#
# La decomposition nx x ny x nz est choisie pour minimiser la surface de
# communication en s'appuyant sur les dimensions physiques du maillage lues
# dans les attributs des balises <lx>, <ly>, <lz> du fichier .arc :
#   <lx nx='150' ...>   <ly ny='150' ...>   <lz nz='150' ...>
# p n'a PAS besoin de diviser exactement le maillage : les sous-domaines
# peuvent etre legerement desequilibres (blocs floor(N/n) ou ceil(N/n)).
#
# Metrique de desequilibre reportee :
#   imbalance = (charge_max - charge_min) / charge_max  in [0, 1[
# Un avertissement est emis si imbalance > seuil (-t, defaut 0.10).
#
# Serie de points de mesure (mode automatique) :
#   diviseurs pairs de Ncores (p impairs exclus), plus p=1.
# Npoints = nombre de ces diviseurs, sauf si -k est fourni.
#
# Usage:
#   ./gen_strong_scaling.sh -i <fichier.arc> -n <Ncores> [options]
#
# Options:
#   -i <fichier>    Fichier .arc de reference (DonneesXXX.arc, balises <nsd>, <lx>, <ly>, <lz>)
#   -n <Ncores>     Nombre maximum de coeurs cible
#   -k <points>     Nombre de points de mesure (defaut: nb diviseurs pairs de Ncores + 1)
#   -o <repertoire> Repertoire de sortie (defaut: ./scaling_strong)
#   -p <liste>      Liste manuelle de valeurs de p separees par des virgules
#                   Ex: -p 1,4,8,16,32,64
#   -t <seuil>      Seuil d'avertissement pour le desequilibre (defaut: 0.10)
#   -v              Mode verbeux
#   -h              Affiche cette aide
#
# Exemples:
#   ./gen_strong_scaling.sh -i DonneesCase1.arc -n 48
#   ./gen_strong_scaling.sh -i DonneesCase1.arc -n 96 -k 6 -t 0.20
#   ./gen_strong_scaling.sh -i DonneesCase1.arc -p 1,4,8,16,32,64
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Valeurs par defaut
# ------------------------------------------------------------------------------
INPUT_FILE=""
NCORES=0
NPOINTS=0          # 0 = automatique
OUTDIR="./scaling_strong"
MANUAL_LIST=""
IMBALANCE_THRESHOLD="0.10"
VERBOSE=0

# ------------------------------------------------------------------------------
# Parsing des arguments
# ------------------------------------------------------------------------------
usage() {
    grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,3\}//' | \
        sed -n '/^Usage/,/^=\{10\}/{ /^=\{10\}/d; p }'
    exit 0
}

while getopts "i:n:k:o:p:t:vh" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG"          ;;
        n) NCORES="$OPTARG"              ;;
        k) NPOINTS="$OPTARG"             ;;
        o) OUTDIR="$OPTARG"              ;;
        p) MANUAL_LIST="$OPTARG"         ;;
        t) IMBALANCE_THRESHOLD="$OPTARG" ;;
        v) VERBOSE=1                     ;;
        h) usage                         ;;
        *) echo "Option invalide: -$OPTARG" >&2; exit 1 ;;
    esac
done

# ------------------------------------------------------------------------------
# Verifications preliminaires
# ------------------------------------------------------------------------------
[[ -z "$INPUT_FILE" ]] && { echo "ERREUR: fichier d'entree requis (-i)" >&2; exit 1; }
[[ ! -f "$INPUT_FILE" ]] && { echo "ERREUR: '$INPUT_FILE' introuvable" >&2; exit 1; }
[[ -z "$MANUAL_LIST" && "$NCORES" -eq 0 ]] && \
    { echo "ERREUR: -n <Ncores> ou -p <liste> requis" >&2; exit 1; }

command -v awk &>/dev/null || { echo "ERREUR: awk est requis" >&2; exit 1; }

log()  { [[ "$VERBOSE" -eq 1 ]] && echo "[INFO] $*" || true; }
warn() { echo "[WARN] $*" >&2; }

# ------------------------------------------------------------------------------
# Extraction des dimensions physiques du maillage depuis les attributs :
#   <lx nx='150' ...>   <ly ny='150' ...>   <lz nz='150' ...>
# ------------------------------------------------------------------------------
MESH_NX=$(grep -oP '<lx[^>]*\snx=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)
MESH_NY=$(grep -oP '<ly[^>]*\sny=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)
MESH_NZ=$(grep -oP '<lz[^>]*\snz=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)

if [[ -z "$MESH_NX" || -z "$MESH_NY" || -z "$MESH_NZ" ]]; then
    echo "ERREUR: impossible de lire nx/ny/nz dans les balises <lx>/<ly>/<lz>" >&2
    echo "       Format attendu : <lx nx='150' ...>  <ly ny='150' ...>  <lz nz='150' ...>" >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Calcul automatique de NPOINTS si non fourni via -k
# Regle : p=1 + diviseurs pairs de Ncores
# En mode liste manuelle, NPOINTS n'est pas utilise.
# ------------------------------------------------------------------------------
if [[ "$NPOINTS" -eq 0 && -z "$MANUAL_LIST" ]]; then
    NPOINTS=$(awk -v n="$NCORES" 'BEGIN{
        count = 1   # p=1 toujours present
        for (d = 2; d <= n; d += 2)
            if (n % d == 0) count++
        print count
    }')
    log "Npoints calcule automatiquement : $NPOINTS (p=1 + diviseurs pairs de $NCORES)"
fi

echo "============================================================"
echo "  Fichier de reference : $INPUT_FILE"
echo "  Ncores               : $NCORES"
echo "  Maillage physique    : $MESH_NX x $MESH_NY x $MESH_NZ"
[[ -z "$MANUAL_LIST" ]] && \
echo "  Points de mesure     : $NPOINTS  (p=1 + diviseurs pairs de $NCORES)"
echo "  Seuil desequilibre   : $IMBALANCE_THRESHOLD"
echo "============================================================"

# ------------------------------------------------------------------------------
# Calcul du desequilibre pour une decomposition (NX/nx, NY/ny, NZ/nz)
#
# Pour chaque axe : les (N mod n) premiers rangs ont ceil(N/n) blocs,
#                   les autres ont floor(N/n) blocs.
# La charge d'un processus est le produit des tailles sur ses 3 axes.
# imbalance = (charge_max - charge_min) / charge_max
# ------------------------------------------------------------------------------
compute_imbalance() {
    local nx=$1 ny=$2 nz=$3   # decomposition processus

    awk -v NX=$MESH_NX -v NY=$MESH_NY -v NZ=$MESH_NZ \
        -v nx=$nx -v ny=$ny -v nz=$nz '
    BEGIN {
        cx_hi = int((NX + nx - 1) / nx)
        cx_lo = int(NX / nx)
        cy_hi = int((NY + ny - 1) / ny)
        cy_lo = int(NY / ny)
        cz_hi = int((NZ + nz - 1) / nz)
        cz_lo = int(NZ / nz)
        load_max = cx_hi * cy_hi * cz_hi
        load_min = cx_lo * cy_lo * cz_lo
        if (load_max == 0) { print "0.0000"; exit }
        printf "%.4f\n", (load_max - load_min) / load_max
    }'
}

# ------------------------------------------------------------------------------
# Surface d'echange totale pour une decomposition.
# Cout = (nx-1)*ny*nz * ceil(NY/ny)*ceil(NZ/nz)
#      + nx*(ny-1)*nz * ceil(NX/nx)*ceil(NZ/nz)
#      + nx*ny*(nz-1) * ceil(NX/nx)*ceil(NY/ny)
# (nombre de faces internes x aire d'une face, approximation conservative)
# ------------------------------------------------------------------------------
surface_cost() {
    local nx=$1 ny=$2 nz=$3
    awk -v NX=$MESH_NX -v NY=$MESH_NY -v NZ=$MESH_NZ \
        -v nx=$nx -v ny=$ny -v nz=$nz '
    BEGIN {
        cx = int((NX + nx - 1) / nx)
        cy = int((NY + ny - 1) / ny)
        cz = int((NZ + nz - 1) / nz)
        cost = (nx-1)*ny*nz * cy*cz \
             + nx*(ny-1)*nz * cx*cz \
             + nx*ny*(nz-1) * cx*cy
        printf "%d\n", cost
    }'
}

# ------------------------------------------------------------------------------
# Meilleure decomposition (nx, ny, nz) pour p processus
# Critere : minimiser surface_cost
# ------------------------------------------------------------------------------
best_decomposition() {
    local p=$1

    local best_nx=1 best_ny=1 best_nz=$p
    local best_cost
    best_cost=$(surface_cost 1 1 $p)

    for (( nx=1; nx<=p; nx++ )); do
        (( p % nx != 0 )) && continue
        local rem=$(( p / nx ))
        for (( ny=1; ny<=rem; ny++ )); do
            (( rem % ny != 0 )) && continue
            local nz=$(( rem / ny ))
            local cost
            cost=$(surface_cost $nx $ny $nz)
            if (( cost < best_cost )); then
                best_cost=$cost
                best_nx=$nx; best_ny=$ny; best_nz=$nz
            fi
        done
    done

    echo "$best_nx $best_ny $best_nz"
}

# ------------------------------------------------------------------------------
# Generation de la serie : {1} union {diviseurs pairs de Ncores}
# Les p impairs (sauf p=1) sont exclus.
# Si -k < nombre de points naturels, sous-echantillonnage geometrique.
# ------------------------------------------------------------------------------
generate_even_divisors_series() {
    local ncores=$1 npts=$2
    local -n _out=$3   # nameref

    local -a natural=(1)
    for (( d=2; d<=ncores; d+=2 )); do
        (( ncores % d == 0 )) && natural+=("$d")
    done
    mapfile -t natural < <(printf '%s\n' "${natural[@]}" | sort -un)
    local total=${#natural[@]}

    if [[ "$npts" -ge "$total" || "$npts" -eq 0 ]]; then
        _out=("${natural[@]}")
    else
        _out=()
        for (( i=0; i<npts; i++ )); do
            local idx
            idx=$(awk -v i=$i -v k=$npts -v n=$total \
                'BEGIN{ printf "%d", int(i*(n-1)/(k-1) + 0.5) }')
            _out+=("${natural[$idx]}")
        done
        mapfile -t _out < <(printf '%s\n' "${_out[@]}" | sort -un)
    fi
}

# ------------------------------------------------------------------------------
# Construction de la liste finale des points p
# ------------------------------------------------------------------------------
declare -a POINTS=()

if [[ -n "$MANUAL_LIST" ]]; then
    IFS=',' read -ra RAW <<< "$MANUAL_LIST"
    for p in "${RAW[@]}"; do
        p=$(echo "$p" | xargs)
        # Filtrage des p impairs (sauf p=1)
        if [[ "$p" -gt 1 && $(( p % 2 )) -ne 0 ]]; then
            warn "p=$p est impair et sera ignore"
            continue
        fi
        POINTS+=("$p")
    done
    mapfile -t POINTS < <(printf '%s\n' "${POINTS[@]}" | sort -un)
else
    generate_even_divisors_series "$NCORES" "$NPOINTS" POINTS
fi

# Garantir p=1
found_one=0
for p in "${POINTS[@]}"; do [[ "$p" -eq 1 ]] && found_one=1; done
[[ "$found_one" -eq 0 ]] && POINTS=(1 "${POINTS[@]}")
mapfile -t POINTS < <(printf '%s\n' "${POINTS[@]}" | sort -un)

# ------------------------------------------------------------------------------
# Repertoire de sortie + fichier recapitulatif
# ------------------------------------------------------------------------------
mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/scaling_plan.txt"

{
echo "# Plan de scalabilite forte"
echo "# Genere le : $(date)"
echo "# Fichier source    : $INPUT_FILE"
echo "# Ncores            : $NCORES"
echo "# Maillage physique : ${MESH_NX}x${MESH_NY}x${MESH_NZ}"
echo "# Seuil desequilibre: $IMBALANCE_THRESHOLD"
echo "#"
printf "# %-6s  %-6s  %-6s  %-6s  %-10s  %s\n" \
    "p" "nx" "ny" "nz" "imbalance" "fichier"
echo "# $(printf '%0.s-' {1..65})"
} > "$SUMMARY"

# ------------------------------------------------------------------------------
# Boucle principale
# ------------------------------------------------------------------------------
echo ""
printf "  %-6s  %-6s  %-6s  %-6s  %-10s  %s\n" \
    "p" "nx" "ny" "nz" "desequil." "fichier"
echo "  $(printf '%0.s-' {1..65})"

GENERATED=0
WARNED=0

for P in "${POINTS[@]}"; do

    # --- Decomposition optimale ---------------------------------------------
    read -r NX NY NZ <<< "$(best_decomposition "$P")"

    # --- Desequilibre -------------------------------------------------------
    IMBALANCE=$(compute_imbalance $NX $NY $NZ)

    OVER_THRESHOLD=$(awk -v im="$IMBALANCE" -v thr="$IMBALANCE_THRESHOLD" \
        'BEGIN{print (im+0 > thr+0) ? 1 : 0}')

    WARN_FLAG=""
    if [[ "$OVER_THRESHOLD" -eq 1 ]]; then
        WARN_FLAG=" !"
        (( WARNED++ )) || true
    fi

    # --- Nom du fichier de sortie -------------------------------------------
    # Format : DonneesXXX_nsd_$(nx)x$(ny)x$(nz).arc  (nx/ny/nz sur 2 chiffres)
    BASE=$(basename "$INPUT_FILE" .arc)
    NX2=$(printf "%02d" "$NX")
    NY2=$(printf "%02d" "$NY")
    NZ2=$(printf "%02d" "$NZ")
    OUTFILE="$OUTDIR/${BASE}_nsd_${NX2}x${NY2}x${NZ2}.arc"

    # --- Generation du fichier .arc -----------------------------------------
    cp "$INPUT_FILE" "$OUTFILE"
    sed -i "s|<nsd>[^<]*</nsd>|<nsd>$NX $NY $NZ</nsd>|" "$OUTFILE"

    printf "  %-6s  %-6s  %-6s  %-6s  %-10s  %s%s\n" \
        "$P" "$NX" "$NY" "$NZ" "$IMBALANCE" "$(basename "$OUTFILE")" "$WARN_FLAG"

    printf "  %-6s  %-6s  %-6s  %-6s  %-10s  %s\n" \
        "$P" "$NX" "$NY" "$NZ" "$IMBALANCE" "$(basename "$OUTFILE")" >> "$SUMMARY"

    (( GENERATED++ )) || true
done

# ------------------------------------------------------------------------------
# Bilan
# ------------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "  Fichiers generes     : $GENERATED"
if [[ "$WARNED" -gt 0 ]]; then
    echo "  !  Points avec desequilibre > $IMBALANCE_THRESHOLD : $WARNED"
    echo "     (ajustez -t pour modifier le seuil d'avertissement)"
fi
echo "  Repertoire de sortie : $OUTDIR/"
echo "  Plan detaille        : $SUMMARY"
echo "============================================================"
