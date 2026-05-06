#!/usr/bin/env bash
# ==============================================================================
# gen_weak_scaling.sh
#
# Genere des fichiers .arc pour une etude de scalabilite faible.
# Pour chaque point de mesure p, le nombre de mailles est multiplie par p
# afin que chaque coeur traite toujours le meme volume de calcul.
#
# Principe :
#   - Maillage de reference : NX0 x NY0 x NZ0 mailles (lu dans <lx>/<ly>/<lz>)
#   - Pour p coeurs, le maillage cible est :
#       NX_p = round(NX0 * p^(1/3))
#       NY_p = round(NY0 * p^(1/3))
#       NZ_p = round(NZ0 * p^(1/3))
#     -> volume total ~ p * NX0*NY0*NZ0, volume par coeur ~ NX0*NY0*NZ0
#   - Le domaine physique est etire isotropiquement par le meme facteur alpha
#     de sorte que la taille d'une maille reste constante :
#       alpha = p^(1/3)
#       contenu de <lx> *= alpha,  contenu de <ly> *= alpha,  contenu de <lz> *= alpha
#   - Les attributs nx/ny/nz de <lx>/<ly>/<lz> sont aussi mis a jour.
#   - La decomposition <nsd>nx ny nz</nsd> est choisie pour minimiser la
#     surface d'echange pour le maillage NX_p x NY_p x NZ_p.
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
#   ./gen_weak_scaling.sh -i <fichier.arc> -n <Ncores> [options]
#
# Options:
#   -i <fichier>    Fichier .arc de reference (DonneesXXX.arc, balises <nsd>, <lx>, <ly>, <lz>)
#   -n <Ncores>     Nombre maximum de coeurs cible
#   -k <points>     Nombre de points de mesure (defaut: nb diviseurs pairs de Ncores + 1)
#   -o <repertoire> Repertoire de sortie (defaut: ./scaling_weak)
#   -p <liste>      Liste manuelle de valeurs de p separees par des virgules
#                   Ex: -p 1,4,8,16,32,64
#   -t <seuil>      Seuil d'avertissement pour le desequilibre (defaut: 0.10)
#   -v              Mode verbeux
#   -h              Affiche cette aide
#
# Exemples:
#   ./gen_weak_scaling.sh -i DonneesCase1.arc -n 48
#   ./gen_weak_scaling.sh -i DonneesCase1.arc -n 96 -k 6 -t 0.20
#   ./gen_weak_scaling.sh -i DonneesCase1.arc -p 1,4,8,16,32,64
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Valeurs par defaut
# ------------------------------------------------------------------------------
INPUT_FILE=""
NCORES=0
NPOINTS=0          # 0 = automatique
OUTDIR="./scaling_weak"
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
# Extraction du maillage de reference depuis les attributs nx/ny/nz
# et des valeurs physiques (contenu texte) de <lx>, <ly>, <lz>
#
# Format attendu (guillemets simples ou doubles) :
#   <lx nx='150' prx='1.0'>1.</lx>
#   <ly ny='150' pry='1.0'>1.</ly>
#   <lz nz='100' prz='1.0'>1.</lz>
# ------------------------------------------------------------------------------
MESH_NX=$(grep -oP '<lx[^>]*\snx=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)
MESH_NY=$(grep -oP '<ly[^>]*\sny=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)
MESH_NZ=$(grep -oP '<lz[^>]*\snz=[\x27"]\K[0-9]+' "$INPUT_FILE" | head -1)

if [[ -z "$MESH_NX" || -z "$MESH_NY" || -z "$MESH_NZ" ]]; then
    echo "ERREUR: impossible de lire nx/ny/nz dans les balises <lx>/<ly>/<lz>" >&2
    echo "       Format attendu : <lx nx='150' ...>  <ly ny='150' ...>  <lz nz='100' ...>" >&2
    exit 1
fi

# Contenu texte : valeur entre > et < (ex: "1." dans <lx nx='150'>1.</lx>)
LX0=$(sed -n 's/.*<lx[^>]*>\([^<]*\)<\/lx>.*/\1/p' "$INPUT_FILE" | head -1 | xargs)
LY0=$(sed -n 's/.*<ly[^>]*>\([^<]*\)<\/ly>.*/\1/p' "$INPUT_FILE" | head -1 | xargs)
LZ0=$(sed -n 's/.*<lz[^>]*>\([^<]*\)<\/lz>.*/\1/p' "$INPUT_FILE" | head -1 | xargs)

if [[ -z "$LX0" || -z "$LY0" || -z "$LZ0" ]]; then
    echo "ERREUR: impossible de lire les valeurs physiques dans <lx>/<ly>/<lz>" >&2
    echo "       Format attendu : <lx nx='150' ...>1.</lx>" >&2
    exit 1
fi

MESH_TOTAL=$(( MESH_NX * MESH_NY * MESH_NZ ))

# ------------------------------------------------------------------------------
# Calcul automatique de NPOINTS
# Regle : p=1 + diviseurs pairs de Ncores
# ------------------------------------------------------------------------------
if [[ "$NPOINTS" -eq 0 && -z "$MANUAL_LIST" ]]; then
    NPOINTS=$(awk -v n="$NCORES" 'BEGIN{
        count = 1
        for (d = 2; d <= n; d += 2)
            if (n % d == 0) count++
        print count
    }')
    log "Npoints calcule automatiquement : $NPOINTS (p=1 + diviseurs pairs de $NCORES)"
fi

echo "============================================================"
echo "  Fichier de reference : $INPUT_FILE"
echo "  Ncores               : $NCORES"
echo "  Maillage reference   : $MESH_NX x $MESH_NY x $MESH_NZ  ($MESH_TOTAL mailles)"
echo "  Domaine physique ref : Lx=$LX0  Ly=$LY0  Lz=$LZ0"
[[ -z "$MANUAL_LIST" ]] && \
echo "  Points de mesure     : $NPOINTS  (p=1 + diviseurs pairs de $NCORES)"
echo "  Seuil desequilibre   : $IMBALANCE_THRESHOLD"
echo "============================================================"

# ------------------------------------------------------------------------------
# Calcul du maillage cible pour p coeurs (etirement isotrope)
#
# alpha = p^(1/3)
# NX_p = round(NX0 * alpha),  NY_p = round(NY0 * alpha),  NZ_p = round(NZ0 * alpha)
# LX_p = LX0 * alpha,         LY_p = LY0 * alpha,         LZ_p = LZ0 * alpha
#
# Sortie : "NX_p NY_p NZ_p LX_p LY_p LZ_p alpha"
# ------------------------------------------------------------------------------
scaled_mesh() {
    local p=$1
    awk -v p=$p \
        -v NX0=$MESH_NX -v NY0=$MESH_NY -v NZ0=$MESH_NZ \
        -v LX0="$LX0"   -v LY0="$LY0"   -v LZ0="$LZ0" '
    BEGIN {
        alpha = p^(1/3)
        NX = int(NX0 * alpha + 0.5)
        NY = int(NY0 * alpha + 0.5)
        NZ = int(NZ0 * alpha + 0.5)
        if (NX < 1) NX = 1
        if (NY < 1) NY = 1
        if (NZ < 1) NZ = 1
        LX = LX0 * alpha
        LY = LY0 * alpha
        LZ = LZ0 * alpha
        printf "%d %d %d %.10g %.10g %.10g %.6f\n", NX, NY, NZ, LX, LY, LZ, alpha
    }'
}

# ------------------------------------------------------------------------------
# Calcul du desequilibre pour une decomposition sur un maillage donne
# imbalance = (charge_max - charge_min) / charge_max
# ------------------------------------------------------------------------------
compute_imbalance() {
    local MNX=$1 MNY=$2 MNZ=$3   # dimensions du maillage
    local nx=$4  ny=$5  nz=$6    # decomposition processus

    awk -v NX=$MNX -v NY=$MNY -v NZ=$MNZ \
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
# Surface d'echange pour un maillage et une decomposition donnes
# ------------------------------------------------------------------------------
surface_cost() {
    local MNX=$1 MNY=$2 MNZ=$3
    local nx=$4  ny=$5  nz=$6
    awk -v NX=$MNX -v NY=$MNY -v NZ=$MNZ \
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
# Meilleure decomposition (nx, ny, nz) pour p processus sur un maillage donne
# Critere : minimiser surface_cost
# ------------------------------------------------------------------------------
best_decomposition() {
    local p=$1 MNX=$2 MNY=$3 MNZ=$4

    local best_nx=1 best_ny=1 best_nz=$p
    local best_cost
    best_cost=$(surface_cost $MNX $MNY $MNZ 1 1 $p)

    for (( nx=1; nx<=p; nx++ )); do
        (( p % nx != 0 )) && continue
        local rem=$(( p / nx ))
        for (( ny=1; ny<=rem; ny++ )); do
            (( rem % ny != 0 )) && continue
            local nz=$(( rem / ny ))
            local cost
            cost=$(surface_cost $MNX $MNY $MNZ $nx $ny $nz)
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
echo "# Plan de scalabilite faible"
echo "# Genere le : $(date)"
echo "# Fichier source    : $INPUT_FILE"
echo "# Ncores            : $NCORES"
echo "# Maillage ref      : ${MESH_NX}x${MESH_NY}x${MESH_NZ}  Lx=$LX0 Ly=$LY0 Lz=$LZ0"
echo "# Seuil desequilibre: $IMBALANCE_THRESHOLD"
echo "#"
printf "# %-4s  %-7s  %-14s  %-14s  %-20s  %-10s  %s\n" \
    "p" "alpha" "maillage" "nsd" "domaine physique" "imbalance" "fichier"
echo "# $(printf '%0.s-' {1..95})"
} > "$SUMMARY"

# ------------------------------------------------------------------------------
# Boucle principale
# ------------------------------------------------------------------------------
echo ""
printf "  %-4s  %-7s  %-14s  %-14s  %-20s  %-10s  %s\n" \
    "p" "alpha" "maillage" "nsd" "domaine physique" "desequil." "fichier"
echo "  $(printf '%0.s-' {1..95})"

GENERATED=0
WARNED=0

for P in "${POINTS[@]}"; do

    # --- Maillage et domaine physique mis a l'echelle -----------------------
    read -r MNX MNY MNZ LX LY LZ ALPHA <<< "$(scaled_mesh "$P")"

    # --- Decomposition optimale pour ce maillage ----------------------------
    read -r NX NY NZ <<< "$(best_decomposition "$P" "$MNX" "$MNY" "$MNZ")"

    # --- Desequilibre -------------------------------------------------------
    IMBALANCE=$(compute_imbalance $MNX $MNY $MNZ $NX $NY $NZ)

    OVER_THRESHOLD=$(awk -v im="$IMBALANCE" -v thr="$IMBALANCE_THRESHOLD" \
        'BEGIN{print (im+0 > thr+0) ? 1 : 0}')

    WARN_FLAG=""
    if [[ "$OVER_THRESHOLD" -eq 1 ]]; then
        WARN_FLAG=" !"
        (( WARNED++ )) || true
    fi

    # --- Nom du fichier de sortie -------------------------------------------
    # Format : DonneesXXX_nsd_NNxNNxNN.arc  (nx/ny/nz sur 2 chiffres)
    BASE=$(basename "$INPUT_FILE" .arc)
    NX2=$(printf "%02d" "$NX")
    NY2=$(printf "%02d" "$NY")
    NZ2=$(printf "%02d" "$NZ")
    OUTFILE="$OUTDIR/${BASE}_nsd_${NX2}x${NY2}x${NZ2}.arc"

    # --- Generation du fichier .arc -----------------------------------------
    cp "$INPUT_FILE" "$OUTFILE"

    # Mise a jour de <nsd>
    sed -i "s|<nsd>[^<]*</nsd>|<nsd>$NX $NY $NZ</nsd>|" "$OUTFILE"

    # Mise a jour des attributs nx/ny/nz dans <lx>/<ly>/<lz>
    sed -i "s|\(<lx[^>]*\)nx='[0-9]*'|\1nx='$MNX'|;
            s|\(<lx[^>]*\)nx=\"[0-9]*\"|\1nx=\"$MNX\"|" "$OUTFILE"
    sed -i "s|\(<ly[^>]*\)ny='[0-9]*'|\1ny='$MNY'|;
            s|\(<ly[^>]*\)ny=\"[0-9]*\"|\1ny=\"$MNY\"|" "$OUTFILE"
    sed -i "s|\(<lz[^>]*\)nz='[0-9]*'|\1nz='$MNZ'|;
            s|\(<lz[^>]*\)nz=\"[0-9]*\"|\1nz=\"$MNZ\"|" "$OUTFILE"

    # Mise a jour du contenu texte de <lx>, <ly>, <lz>
    sed -i "s|<lx\([^>]*\)>[^<]*</lx>|<lx\1>$LX</lx>|" "$OUTFILE"
    sed -i "s|<ly\([^>]*\)>[^<]*</ly>|<ly\1>$LY</ly>|" "$OUTFILE"
    sed -i "s|<lz\([^>]*\)>[^<]*</lz>|<lz\1>$LZ</lz>|" "$OUTFILE"

    MESH_STR="${MNX}x${MNY}x${MNZ}"
    NSD_STR="${NX}x${NY}x${NZ}"
    DOM_STR="$(printf "%.4g" $LX)x$(printf "%.4g" $LY)x$(printf "%.4g" $LZ)"

    printf "  %-4s  %-7s  %-14s  %-14s  %-20s  %-10s  %s%s\n" \
        "$P" "$ALPHA" "$MESH_STR" "$NSD_STR" "$DOM_STR" \
        "$IMBALANCE" "$(basename "$OUTFILE")" "$WARN_FLAG"

    printf "  %-4s  %-7s  %-14s  %-14s  %-20s  %-10s  %s\n" \
        "$P" "$ALPHA" "$MESH_STR" "$NSD_STR" "$DOM_STR" \
        "$IMBALANCE" "$(basename "$OUTFILE")" >> "$SUMMARY"

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
