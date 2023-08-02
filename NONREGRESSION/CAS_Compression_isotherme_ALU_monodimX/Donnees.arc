<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Exemple Arcane d'un module Hydro tres simplifie</title>
    <timeloop>MahycoLoop</timeloop>
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
    <output-history-period>500</output-history-period>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>1 1</nsd> 
       <origine>-5.e-6 -5.e-6</origine>
       <lx nx='100' prx='1.0'>10.e-6</lx>
       <ly ny='10' pry='1.0'>10.e-6</ly>
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
  <material><name>Al_mat</name></material>
  <environment>
    <name>ALU</name>
    <material>Al_mat</material>
    <densite-initiale>2785.00</densite-initiale>
    <pression-initiale>0.</pression-initiale>
    <eos-model name="MieGruneisen">
      <adiabatic-cst>2.</adiabatic-cst>
      <rho0>2785.</rho0>
      <c-cst>7.905e10</c-cst>
      <d-cst>1.3266e11</d-cst>
      <s-cst>0.</s-cst>
      <specific-heat>1000.</specific-heat>
    </eos-model>
    <linear-pseudo-coeff>0.5</linear-pseudo-coeff>
    <quadratic-pseudo-coeff>1.2</quadratic-pseudo-coeff>
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
     <deltat-init>1.e-12</deltat-init>
     <deltat-min>1.e-14</deltat-min>
     <deltat-max>1.e-12</deltat-max>
    <longueur-caracteristique>monodimX</longueur-caracteristique>
    <cfl>0.1</cfl>
    <final-time>1.e-7</final-time>
    
    <time-history>
    <periode>10</periode>
    <borne-inf>0.e-6 0.e-6 -1.</borne-inf>
    <borne-sup>2.e-6 2.e-6 1.</borne-sup>
    </time-history>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>LinearPressure</type>
      <value>0.</value>
      <dependanceT>0.3e18</dependanceT>
      <Tdebut>0.e-9</Tdebut>   
      <Tfin>100.e-9</Tfin>    
    </boundary-condition>

    <boundary-condition>
      <surface>XMIN</surface>
      <type>GeometricPressure</type>
      <value>30.</value>
      <Tdebut>100.e-9</Tdebut>   
      <Tfin>200.e-9</Tfin>       
    </boundary-condition>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>Vy</type>
      <value>0.</value>
    </boundary-condition>
    
    <boundary-condition>
      <surface>XMAX</surface>
      <type>LinearPressure</type>
      <value>0.</value>
      <dependanceT>0.3e18</dependanceT>
      <Tdebut>0.e-9</Tdebut>   
      <Tfin>100.e-9</Tfin>         
    </boundary-condition>

    <boundary-condition>
      <surface>XMAX</surface>
      <type>GeometricPressure</type>
      <value>30e9</value>
      <Tdebut>100.e-9</Tdebut>   
      <Tfin>200.e-9</Tfin>      
    </boundary-condition>
    
    <boundary-condition>
      <surface>XMAX</surface>
      <type>Vy</type>
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
