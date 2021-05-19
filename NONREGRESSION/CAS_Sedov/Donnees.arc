<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Exemple Arcane d'un module Hydro tres simplifie</title>
    <timeloop>MahycoLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>100</output-period>
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
       <nsd>2 2 2</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='60' prx='1.0'>1.2</lx>

       <ly ny='60' pry='1.0'>1.2</ly>

       <lz nz='60' prz='1.0'>1.2</lz>
     </cartesian>

     </meshgenerator>

    <initialisation>
    </initialisation>
  </mesh>

  <arcane-checkpoint>
    <period>500</period>
    <!-- Mettre '0' si on souhaite ne pas faire de protections a la fin du calcul -->
    <do-dump-at-end>0</do-dump-at-end>
    <checkpoint-service name="ArcaneBasic2CheckpointWriter" />
  </arcane-checkpoint>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  <material><name>Init</name></material>
  <material><name>Air</name></material>
  <environment>
    <name>ZoneInit</name>
    <material>Init</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
    </eos-model> 
  </environment>
  <environment>
    <name>ZoneAir</name>
    <material>Air</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
    <!-- <eos-model name="StiffenedGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <limit-tension>0.01</limit-tension> -->
    </eos-model> 
  </environment>
   
   <cas-test>11</cas-test>
    <pseudo-centree>0</pseudo-centree>
    <ordre-projection>2</ordre-projection>
    <schema-csts>0</schema-csts>
     <deltat-init>0.00001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.01</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
    <convnerg>false</convnerg>
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
    <boundary-condition>
      <surface>ZMIN</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
  </mahyco>
</case>
