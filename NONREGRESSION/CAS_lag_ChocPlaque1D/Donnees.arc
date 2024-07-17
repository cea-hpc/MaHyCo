<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_lag_ChocPlaque1D</title>
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
       <lx nx='100' prx='1.0'>2.e-6</lx>
       <ly ny='10' pry='1.0'>1.e-7</ly>
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
  </time-history>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
  <material><name>AL_mat</name></material>
  <environment>
    <name>ALU</name>
    <material>AL_mat</material>
    <!--<densite-initiale>7286.62</densite-initiale>-->
    <densite-initiale>7290.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
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
     <deltat-init>5.e-15</deltat-init>
     <deltat-min>1.e-16</deltat-min>
     <deltat-max>5.e-15</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>1.e-12</final-time>
    
    <time-history>
    <periode>5000</periode>
    <borne-inf>0.8e-6 -1 -1.</borne-inf>
    <borne-sup>0.82e-6 1. 1.</borne-sup>
    <borne-inf-noeud>1.99e-6 -1 -1.</borne-inf-noeud>
    <borne-sup-noeud>2.01e-6 1. 1.</borne-sup-noeud>
    </time-history>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vx</type>
      <value>3600</value>
      <Tdebut>0.e-9</Tdebut>   
      <Tfin>0.02e-9</Tfin>         
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
