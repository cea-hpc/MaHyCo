<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_2D_Sedov</title>
    <timeloop>MahycoLoop</timeloop>
    <modules>
      <module name="TimeHistory" active="false" />
    </modules>
  </arcane>

  <arcane-post-processing>
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
       <nsd>1 1</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='60' prx='1.0'>1.2</lx>
       <ly ny='60' pry='1.0'>1.2</ly>
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
  <!-- <bilan name="NodeWatching">
  </bilan> -->
  <!-- <bilan name="CellWatching">
  </bilan> -->
    <bilan name="EnvSummation">
    <!-- <environment>ZoneAir</environment> -->
  </bilan>
  </time-history>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  <material><name>Air</name></material>
  <material><name>Init</name></material>
  <environment>
    <name>ZoneAir</name>
    <material>Air</material>
    <densite-initiale>1.</densite-initiale>
    <pression-initiale>0.0979264e-5</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
    <!-- <eos-model name="StiffenedGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
      <limit-tension>0.01</limit-tension> -->
    </eos-model> 
  </environment>
  <environment>
    <name>ZoneInit</name>
    <material>Init</material>
    <densite-initiale>1.</densite-initiale>
    <pression-initiale>1.</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
    </eos-model> 
  </environment>

  
   <cas-model name="SEDOV">
   <cas-test>1</cas-test>
   </cas-model>
   <remap name="RemapADI">
    <!-- <is-euler-scheme>true</is-euler-scheme>
    <volum-criteria>0.8</volum-criteria>
    <nb-iteration-winslow>10</nb-iteration-winslow> -->
    <ordre-projection>2</ordre-projection>
   </remap>
       <!--<with-projection>false</with-projection>-->
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.00001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.01</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>1.</final-time>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
  </mahyco>
</case>
