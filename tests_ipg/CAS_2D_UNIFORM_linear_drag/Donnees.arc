<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_2D_UNIFORM</title>
    <timeloop>MahycoLoop</timeloop>
    <!-- <timeloop>MahycoLagrangeLoop</timeloop> -->
    <modules>
      <module name="ArcanePostProcessing" active="true" />
      <module name="Ipg" active="true"/>
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
    <output-frequency>1e5</output-frequency>
    <format>
      <binary-file>true</binary-file>
    </format>
  </arcane-post-processing>
  
  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>1 1</nsd>
       <origine>0. 0. 0.</origine>
       <lx nx='10' prx='1.0'>1e5</lx>
       <ly ny='1' pry='1.0'>1e4</ly>
     </cartesian>
     </meshgenerator>
    <initialisation>
    </initialisation>
  </mesh>

  <arcane-checkpoint>
    <period>1000</period>
    <!-- Mettre '0' si on souhaite ne pas faire de protections a la fin du calcul -->
    <do-dump-at-end>0</do-dump-at-end>
    <checkpoint-service name="ArcaneBasic2CheckpointWriter" />
  </arcane-checkpoint>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>

    <material><name>gaz_mat</name></material>
    <environment>
      <name>gaz</name>
      <material>gaz_mat</material>
      <densite-initiale>1.</densite-initiale>
      <pression-initiale>250.</pression-initiale>
      <vitesse-initiale>0. 0. 0.</vitesse-initiale>
      <temperature-initiale>1.20267</temperature-initiale>
      <eos-model name="PerfectGas">
        <adiabatic-cst>1.4</adiabatic-cst>
        <specific-heat>1.</specific-heat>
      </eos-model> 
    </environment>
    
    <cas-model name="UNIFORM">
    </cas-model>

    <pseudo-centree>0</pseudo-centree>
    <with-projection>false</with-projection>
    <schema-csts>0</schema-csts>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>

    <!-- temps final et timestep -->
    <deltat-init>1e-1</deltat-init>
    <deltat-min>1e-3</deltat-min>
    <deltat-max>1e4</deltat-max>
    <!-- <cfl>0.05</cfl> -->
    <final-time>1e6</final-time>

    <!-- BC -->
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
  

  <!-- IPG -->
  <ipg>
    <particle-injector-type name="OneParticle">
      <init-velocity>1. 0. 0.</init-velocity>
      <init-coord>.05 4000. 0.</init-coord>
      <init-time>0.</init-time>
      <init-radius>1.</init-radius>
      <init-temperature>1.</init-temperature>
      <init-density>4125.</init-density>
    </particle-injector-type>
    <ipg-output name="VtkAscii">
      <variable>ParticleVelocity</variable>
      <variable>ParticleRadius</variable>
      <variable>ParticleTemperature</variable>
      <output-frequency>1e5</output-frequency>
    </ipg-output>
    <spray-type name="SprayFin">
      <drag>
        <type>Linear</type>
        <coef>3.85e-5</coef>
      </drag>
    </spray-type>
    <one-particle-history>1</one-particle-history>
  </ipg>

  
</case>
