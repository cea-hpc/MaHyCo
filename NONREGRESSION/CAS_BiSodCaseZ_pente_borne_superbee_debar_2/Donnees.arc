<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Exemple Arcane d'un module Hydro tres simplifie</title>
    <timeloop>MahycoLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>50</output-period>
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
       <nsd>1 1 4</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='5' prx='1.0'>.1</lx>

       <ly ny='5' pry='1.0'>.1</ly>

       <lz nz='100' prz='1.0'>1.</lz>
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
  <material><name>ZG_mat</name></material>
  <material><name>ZD_mat</name></material>
  <environment>
    <name>ZG</name>
    <material>ZG_mat</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
    </eos-model> 
  </environment>
  <environment>
    <name>ZD</name>
    <material>ZD_mat</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
    <!-- <eos-model name="StiffenedGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <limit-tension>0.01</limit-tension> -->
    </eos-model> 
  </environment>
   
   <cas-model name="SOD">
   <cas-test>15</cas-test>
   </cas-model>
    <pseudo-centree>0</pseudo-centree>
    <ordre-projection>2</ordre-projection>
    <projection-pente-borne>true</projection-pente-borne>
    <projection-pente-borne-debar-fix>2</projection-pente-borne-debar-fix>
    <projection-limiteur-id>1</projection-limiteur-id>
    <schema-csts>0</schema-csts>
     <deltat-init>0.00001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.01</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
    <convnerg>false</convnerg>
    <final-time>.2</final-time>
    
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
    <boundary-condition>
      <surface>ZMIN</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>ZMAX</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
		
  </mahyco>
</case>
