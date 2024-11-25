<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_shallow</title>
    <timeloop>MahycoLoop</timeloop>
    <modules>
      <module name="TimeHistory" active="true" />
    </modules>
  </arcane>

  <arcane-post-processing>
    <save-init>true</save-init>
    <output>
      <variable>CellMass</variable>
      <variable>Pressure</variable>
      <variable>Density</variable>
      <variable>Velocity</variable>
      <variable>NodeMass</variable>
      <variable>InternalEnergy</variable>
      <variable>PseudoViscosity</variable>
      <variable>Materiau</variable>
      <variable>Temperature</variable>
      <variable>SoundSpeed</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
  </arcane-post-processing>

  <mesh nb-ghostlayer="2" ghostlayer-builder-version="3">
    <meshgenerator>
     <cartesian>
       <nsd>4 2</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='100' prx='1.0'>1.</lx>

       <ly ny='50' pry='1.0'>.5</ly>
     </cartesian>

     </meshgenerator>

    <initialisation>
    </initialisation>
  </mesh>

  <arcane-checkpoint>
    <period>4000</period>
    <!-- Mettre '0' si on souhaite ne pas faire de protections a la fin du calcul -->
    <do-dump-at-end>0</do-dump-at-end>
    <checkpoint-service name="ArcaneBasic2CheckpointWriter" />
  </arcane-checkpoint>

  <time-history>
    <bilan name="EnvSummation">
      <variable>CellMass</variable>
      <variable>Pressure</variable>
      <variable>Density</variable>
      <variable>InternalEnergy</variable>
    </bilan>
  </time-history>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  
  <material><name>Fictif_mat</name></material>
  <material><name>Plaque_mat</name></material>
  <material><name>Bulle_mat</name></material>
  <material><name>Fictif2_mat</name></material>
  
  <environment>
    <name>Fictif</name>
    <material>Fictif_mat</material>
    <densite-initiale>2785.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <eos-model name="Fictif">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>1004</specific-heat>
      <tdebut-pression>0.</tdebut-pression>
      <tfin-pression>.25</tfin-pression>
      <valeur-debut-pression>3.e9</valeur-debut-pression>
      <valeur-dependance-temps>0.</valeur-dependance-temps>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
  </environment>
  
  <environment>
    <name>Plaque</name>
    <material>Plaque_mat</material>
    <densite-initiale>2785.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <eos-model name="MieGruneisen">
     <adiabatic-cst>2.</adiabatic-cst>
      <rho0>2785.</rho0>
      <s1>1.339</s1>
      <s2>0.</s2>
      <s3>0.</s3>
      <gamma0>1.6</gamma0>
      <a>0.48</a>
      <c0>5328.</c0>
      <energie-ref>0.</energie-ref>
      <specific-heat>897.</specific-heat>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
    <elasto-model name="NoModel">
    </elasto-model>
    <linear-pseudo-coeff>0.5</linear-pseudo-coeff>
    <quadratic-pseudo-coeff>1.</quadratic-pseudo-coeff>
  </environment>
  
  <environment>
    <name>Bulle</name>
    <material>Bulle_mat</material>
    <densite-initiale>.001</densite-initiale>
    <pression-initiale>1.e2</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>1004</specific-heat>
    </eos-model> 
  </environment>
  
  <environment>
    <name>Fictif2</name>
    <material>Fictif2_mat</material>
    <densite-initiale>.001</densite-initiale>
    <pression-initiale>1.e2</pression-initiale>
    <eos-model name="Fictif">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>1004</specific-heat>
      <tdebut-pression>0.</tdebut-pression>
      <tfin-pression>100000.</tfin-pression>
      <valeur-debut-pression>1.e5</valeur-debut-pression>
      <valeur-dependance-temps>0.</valeur-dependance-temps>
    </eos-model>
  </environment>
  
   <cas-model name="CHOCBULLE">
   <cas-test>47</cas-test>
   <nombre-bulles-x>2</nombre-bulles-x>
   <nombre-bulles-y>2</nombre-bulles-y>
   <position-premiere-bulle>0.3 0.1 0.0</position-premiere-bulle>
   <deltax-bulle>0.1</deltax-bulle>
   <deltay-bulle>0.1</deltay-bulle>
   <rayon-bulle>0.02</rayon-bulle>
   <nombre-bulles-aleatoires>0</nombre-bulles-aleatoires>
   <position-min>0.35 0.25 0.</position-min>
   <position-max>0.5 0.45 0.</position-max>
   <rayon-min>0.01</rayon-min>
   <rayon-max>0.03</rayon-max>
   </cas-model>
   <remap name="RemapADI">
    <ordre-projection>2</ordre-projection>
   </remap>
   
    <pseudo-centree>0</pseudo-centree>
    <with-newton>false</with-newton>
    <pression-explicite>true</pression-explicite>
    <schema-csts>0</schema-csts>
     <deltat-init>1.e-8</deltat-init>
     <deltat-min>1.e-11</deltat-min>
     <deltat-max>0.01</deltat-max>
     <cfl>0.05</cfl>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>.00005</final-time>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>XMAX</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMAX</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
		
  </mahyco>
</case>
