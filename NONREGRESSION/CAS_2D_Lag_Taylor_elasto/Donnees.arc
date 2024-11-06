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
      <variable>SoundSpeed</variable>
      <variable>VelocityGradient</variable>
      <variable>DeformationRate</variable>
      <variable>StrainTensorXX</variable>
      <variable>StrainTensorXY</variable>
      <variable>StrainTensorYY</variable>
      <variable>PlasticDeformationVelocity</variable>
      <variable>PlasticDeformation</variable>

    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
  </arcane-post-processing>

  <!-- <meshes> -->
	<mesh>
        <file internal-partition='true'>plaqueTaylor.msh</file>
      <initialisation>
        <variable nom="Density" valeur="8129." groupe="Metal"/>
        <variable nom="Pressure" valeur="1.e5" groupe="Metal"/>
        <variable nom="Materiau" valeur="0." groupe="Metal" />
      </initialisation>
    </mesh>
    <!-- </meshes> -->

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
  <material><name>Metal</name></material>
  <environment>
    <name>Acier</name>
    <material>Metal</material>
    <vitesse-initiale>0. -144. 0.</vitesse-initiale>
    <initialisation-utilisateur>true</initialisation-utilisateur>
    <eos-model name="MieGruneisen">
      <adiabatic-cst>1.6</adiabatic-cst>
      <rho0>8129.</rho0>
      <s1>1.58</s1>
      <s2>0.</s2>
      <s3>0.</s3>
      <gamma0>1.6</gamma0>
      <a>0.5</a>
      <c0>3980.</c0>
      <energie-ref>0.</energie-ref>
      <specific-heat>502.</specific-heat>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
    <linear-pseudo-coeff>0.5</linear-pseudo-coeff>
    <quadratic-pseudo-coeff>1.2</quadratic-pseudo-coeff>
    <elasto-model name="DefaultModel">
      <elasto-mu-model name="EPP">
        <elastic-cst>6.4e10</elastic-cst>
      </elasto-mu-model>
      <elasto-y-model name="YEPP">
        <limit-elastic>3.e8</limit-elastic>
      </elasto-y-model>
    </elasto-model> 
  </environment>
  
    <with-projection>false</with-projection>
    <coef-antiderive>0.01</coef-antiderive>
    <pseudo-centree>0</pseudo-centree>
    <schema-csts>0</schema-csts>
     <deltat-init>0.0000001</deltat-init>
     <deltat-min>0.0000000001</deltat-min>
     <deltat-max>0.0001</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
    <cfl>0.05</cfl> 
    <final-time>3.e-5</final-time>
    
   
    <boundary-condition>
      <surface>YMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>		
  </mahyco>
</case>
