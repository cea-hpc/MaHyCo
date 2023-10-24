<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Noh</title>
    <timeloop>MahycoLoop</timeloop>
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

  <mesh>
      <file internal-partition='true'>noh.msh</file>
      <!-- <file internal-partition='true'>noh.unv</file>-->
      <!-- <file internal-partition='true'>noh.vtk</file>-->
      <initialisation>
        <variable nom="Density" valeur="1." groupe="AIR"/>
        <variable nom="Pressure" valeur="1.e-14" groupe="AIR"/>
        <variable nom="Materiau" valeur="0." groupe="AIR" />
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
  <material><name>AIR</name></material>
  
  <environment>
    <name>NOH</name>
    <material>AIR</material>
    <initialisation-utilisateur>true</initialisation-utilisateur>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.66</adiabatic-cst>
      <specific-heat>12.47</specific-heat>
    </eos-model> 
  </environment>
 
    <with-projection>false</with-projection>
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.00001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.0001</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>.6</final-time>
    
    <boundary-condition>
      <surface>RMIN</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>RMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>RMAX</surface>
      <type>OFV</type>
      <value>1.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
  </mahyco>
</case>
