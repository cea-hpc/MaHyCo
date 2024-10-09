<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Sod</title>
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

  <!-- <meshes> -->
	<mesh>
        <!-- <filename>sod.unv</filename>-->
        <!-- <file>sod.msh</file>-->
        <file internal-partition='true'>sod.msh</file>
      <initialisation>
        <variable nom="Density" valeur="1." groupe="ZG_mat"/>
        <variable nom="Pressure" valeur="1." groupe="ZG_mat"/>
        <variable nom="Materiau" valeur="0." groupe="ZG_mat" />
        <variable nom="Density" valeur="0.125" groupe="ZD_mat"/>
        <variable nom="Pressure" valeur=".1" groupe="ZD_mat"/>
        <variable nom="Materiau" valeur="1." groupe="ZD_mat" />
      </initialisation>
    </mesh>
    <!-- nsd 4 -->

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
  <material><name>ZG_mat</name></material>
  <material><name>ZD_mat</name></material>
  <environment>
    <name>ZG</name>
    <material>ZG_mat</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>20.785</specific-heat>
    </eos-model> 
  </environment>
  <environment>
    <name>ZD</name>
    <material>ZD_mat</material>
    <eos-model name="PerfectGas">
      <adiabatic-cst>1.4</adiabatic-cst>
      <specific-heat>20.785</specific-heat>
    </eos-model> 
  </environment>
  
    <with-projection>false</with-projection>
    <coef-antiderive>0.01</coef-antiderive>
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.00001</deltat-init>
     <deltat-min>0.00000001</deltat-min>
     <deltat-max>0.0001</deltat-max>
    <longueur-caracteristique>monodimX</longueur-caracteristique>
     
    <final-time>.2</final-time>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>XMAX</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMIN</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>YMAX</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
		
  </mahyco>
</case>
