<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Exemple Arcane d'un module Hydro tres simplifie</title>
    <timeloop>MahycoLoop</timeloop>
  </arcane>

  <arcane-post-processing>
  <save-init>true</save-init>
    <output-period>50</output-period>
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
      <variable>PlasticDeformationVelocity</variable>
      <variable>PlasticDeformation</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
    <output-history-period>500</output-history-period>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>1 1</nsd> 
       <origine>0.0 0.0</origine>
       <lx nx='100' prx='1.0'>100.e-6</lx>

       <ly ny='50' pry='1.0'>50.e-6</ly>
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
  <material><name>AL_mat</name></material>
  <environment>
    <name>ALU</name>
    <material>AL_mat</material>
    <!--<densite-initiale>7286.62</densite-initiale>-->
    <densite-initiale>7290.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <!--<energie-initiale>7.738376</energie-initiale>-->
    <!--<energie-initiale>7.47539332</energie-initiale>-->
    <!-- <energie-initiale>7.8112</energie-initiale> donne 1e5 en pression -->
    <!-- <energie-initiale>-1.7090e+3</energie-initiale>  celui du lagrange à 1e5 avec la sesame -->
    <energie-initiale>6.091459</energie-initiale>  <!-- celui de l'euler à 2.6593 avec le lagrange en sesame donnant P=2.65e+7 -->
    <eos-model name="Autre">
        <fichier-coeff>ee.CineTest22#.Sn.00#.coeff</fichier-coeff>
    </eos-model> 
    <elasto-model name="NoModel">
        <elastic-cst>3.e10</elastic-cst>
        <limit-elastic>1.e9</limit-elastic>
    </elasto-model> 
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
     
    <final-time>1.e-7</final-time>
    <!--<final-time>2.e-10</final-time>-->
    
    <time-history>
    <periode>50</periode>
    <borne-sup>91.e-6 6.e-6 0.</borne-sup>
    <borne-inf>90.e-6 5.e-6 0.</borne-inf>
    </time-history>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>GeometricPressure</type>
      <value>10.e9</value>
      <dependanceY>2.e12</dependanceY> <!-- Y varie de 0.e-6 à 50e-6 donc pression de 10e9 à 20e9 -->
      <Tdebut>1.e-11</Tdebut>   
      <Tfin>9.01e-9</Tfin>   
    </boundary-condition>
    <boundary-condition>
      <surface>XMIN</surface>
      <type>GeometricPressure</type>
      <value>1.e5</value>
      <Tdebut>10.e-9</Tdebut>   
      <Tfin>100.e-9</Tfin>      
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
