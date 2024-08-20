<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_MonoRiderDeformationTimeReverse</title>
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
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
  </arcane-post-processing>

  <mesh nb-ghostlayer="2" ghostlayer-builder-version="3">
 
    <meshgenerator>
     <cartesian>
       <nsd>4 2 1</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='50' prx='1.0'>1.0</lx>

       <ly ny='50' pry='1.0'>1.0</ly>

       <lz nz='1' prz='1.0'>0.1</lz>
     </cartesian>
     </meshgenerator>
    <initialisation>
    </initialisation>
  </mesh>

  <arcane-checkpoint>
    <period>0</period>
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
  <material><name>Vide_mat</name></material>
  <material><name>Bulle_mat</name></material>
  <environment>
    <name>Vide</name>
    <material>Vide_mat</material>
    <densite-initiale>0.</densite-initiale>
    <pression-initiale>0.</pression-initiale>
    <eos-model name="Fictif">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
      <tdebut-pression>0.</tdebut-pression>
      <tfin-pression>10</tfin-pression>
      <valeur-debut-pression>0.</valeur-debut-pression>
      <valeur-dependance-temps>0.</valeur-dependance-temps>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
  </environment>
  <environment>
    <name>Bulle</name>
    <material>Bulle_mat</material>
    <densite-initiale>1.</densite-initiale>
    <pression-initiale>0.</pression-initiale>
    <eos-model name="Fictif">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
      <tdebut-pression>0.</tdebut-pression>
      <tfin-pression>10</tfin-pression>
      <valeur-debut-pression>0.</valeur-debut-pression>
      <valeur-dependance-temps>0.</valeur-dependance-temps>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
  </environment>
   
   <cas-model name="RIDER">
   <cas-test>30</cas-test>
   <reverse-option>true</reverse-option>
   <parametre>1.</parametre>
   </cas-model>
   <remap name="RemapADI">
    <ordre-projection>2</ordre-projection>
   </remap>
   
    <sans-lagrange>true</sans-lagrange>
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.0001</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>1.</final-time>
    
 <!--   <boundary-condition>
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
    <boundary-condition>
      <surface>ZMIN</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>ZMAX</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>-->
		
  </mahyco>
</case>
