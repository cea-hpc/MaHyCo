<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_MonoPointtriple_pb_sp_mixte_superbeeG_superbeeG</title>
    <timeloop>MahycoLoop</timeloop>
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
       <nsd>2 2</nsd> 
       <origine>0.0 0.0</origine>
       <lx nx='70' prx='1.0'>.07</lx>

       <ly ny='30' pry='1.0'>.03</ly>
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

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  <material><name>Mat1</name></material>
  <material><name>Mat2</name></material>
  <material><name>Mat3</name></material>
  <environment>
    <name>Mat1</name>
    <material>Mat1</material>
    <densite-initiale>1.</densite-initiale>
    <pression-initiale>1.</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
    </eos-model> 
  </environment>
  <environment>
    <name>Mat2</name>
    <material>Mat2</material>
    <densite-initiale>1.</densite-initiale>
    <pression-initiale>0.1</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
    </eos-model> 
  </environment>
   <environment>
    <name>Mat3</name>
    <material>Mat3</material>
    <densite-initiale>.1</densite-initiale>
    <pression-initiale>0.1</pression-initiale>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>2.4</specific-heat>
    </eos-model> 
  </environment>
  
   <cas-model name="POINTTRIPLE">
   <cas-test>12</cas-test>
   </cas-model>
   <with-projection>true</with-projection>
   <remap name="RemapADI">
    <ordre-projection>2</ordre-projection>
    <projection-pente-borne>true</projection-pente-borne>
    <projection-simple-pente>true</projection-simple-pente>
    
    <projection-pente-borne-mixte>true</projection-pente-borne-mixte>
    <projection-pente-borne-debar-fix>1</projection-pente-borne-debar-fix>
    <projection-limiteur-id>4</projection-limiteur-id>
    <projection-limiteur-pure-id>4</projection-limiteur-pure-id>
    <threshold>1.e-16</threshold>
   </remap>
    <sans-lagrange>false</sans-lagrange>
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.000000001</deltat-init>
     <deltat-min>0.0000000001</deltat-min>
     <deltat-max>0.00001</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
    <threshold>1.e-16</threshold>
    <final-time>.041</final-time>
    
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
