<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_lag_ChocPlaque3D</title>
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
      <variable>Temperature</variable>
      <variable>Density</variable>
      <variable>Velocity</variable>
      <variable>NodeMass</variable>
      <variable>InternalEnergy</variable>
      <variable>PseudoViscosity</variable>
      <variable>Materiau</variable>
      <variable>FracPhase1</variable>
      <variable>FracPhase2</variable>
      <variable>FracPhase3</variable>
      <variable>VelocityGradient</variable>
      <variable>DeformationRate</variable>
      <variable>StrainTensorXX</variable>
      <variable>StrainTensorXY</variable>
      <variable>StrainTensorXZ</variable>
      <variable>StrainTensorYY</variable>
      <variable>StrainTensorYZ</variable>
      <variable>PlasticDeformationVelocity</variable>
      <variable>PlasticDeformation</variable>
      <variable>StrainNorm</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>2 2 2</nsd> 
       <origine>-20.e-6 -20.e-6 0.</origine>
       <lx nx='60' prx='1.0'>40.e-6</lx>
       <ly ny='60' pry='1.0'>40.e-6</ly>
       <lz nz='30' prz='1.0'>20.e-6</lz>
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

  <time-history>
    <bilan name="EnvSummation">
      <variable>CellMass</variable>
      <variable>Pressure</variable>
      <variable>Density</variable>
      <variable>InternalEnergy</variable>
    </bilan>
     <bilan name="CellWatching">
     <borne-inf>0.e-6 0.e-6 1.8e-5</borne-inf>
     <borne-sup>0.5e-6 0.5e-6 1.84e-5</borne-sup>
    </bilan>
  </time-history>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  <material><name>AL_mat</name></material>
  <environment>
    <name>ALU</name>
    <material>AL_mat</material>
    <densite-initiale>7290.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <energie-initiale>6.091459</energie-initiale>  <!-- celui de l'euler à 2.6593 avec le lagrange en sesame donnant P=2.65e+7 -->
    <eos-model name="Autre">
        <fichier-coeff>ee.CineTest22#.Sn.00#.coeff</fichier-coeff>
    </eos-model> 
    <elasto-model name="DefaultModel">
     <yandg-model name="EPP">
        <elastic-cst>3.e10</elastic-cst>
        <limit-elastic>1.e9</limit-elastic>
     </yandg-model>
    </elasto-model>
    <linear-pseudo-coeff>0.15 </linear-pseudo-coeff>
    <quadratic-pseudo-coeff>0.5 </quadratic-pseudo-coeff>
  </environment>
  
   
   <cas-model name="ONECELL">
   <cas-test>0</cas-test>
   </cas-model>
   
    <pseudo-centree>0</pseudo-centree>
    <with-newton>false</with-newton>
    <with-projection>false</with-projection>
    <pression-explicite>true</pression-explicite>
    <schema-csts>0</schema-csts>
     <deltat-init>1.e-11</deltat-init>
     <deltat-min>1.e-14</deltat-min>
     <deltat-max>1.e-11</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>5.e-10</final-time>
    
    <boundary-condition>
      <surface>ZMAX</surface>
      <type>GeometricPressure</type>
      <value>10.e9</value>
      <Tdebut>1.e-11</Tdebut>   
      <Tfin>49.01e-9</Tfin>   
    </boundary-condition>
    <boundary-condition>
      <surface>ZMIN</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>XMAX</surface>
      <type>Vx</type>
      <value>0.</value>
    </boundary-condition>
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
      <surface>YMAX</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
		
  </mahyco>
</case>
