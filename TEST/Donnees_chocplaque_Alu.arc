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
       <lx nx='500' prx='1.0'>500.e-6</lx>

       <ly ny='10' pry='1.0'>10.e-6</ly>
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
  <material><name>AL_mat</name></material>
  <environment>
    <name>ALU</name>
    <material>AL_mat</material>
    <densite-initiale>2710.</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <eos-model name="MieGruneisen">
 <!--   <adiabatic-cst>2.1</adiabatic-cst>
        <specific-heat>1000.</specific-heat>
        <rho0>2710.</rho0>
        <c-cst>7.84e+10</c-cst>
        <d-cst>4.89e+10</d-cst>
        <s-cst>-6.04e+10</s-cst>--> 
        
        <adiabatic-cst>2.</adiabatic-cst>
        <specific-heat>1000.</specific-heat>
        <rho0>2710.</rho0>
        <c-cst>7.905e+10</c-cst>
        <d-cst>1.3266e+11</d-cst>
        <s-cst>0.0</s-cst>
    </eos-model> 
    <elasto-model name="DefaultModel">
        <elastic-cst>0.e10</elastic-cst>
        <limit-elastic>1.e9</limit-elastic>
    </elasto-model> 
  </environment>
  
   
   <cas-model name="ONECELL">
   <cas-test>0</cas-test>
   </cas-model>
   
    <pseudo-centree>0</pseudo-centree>
    <with-newton>true</with-newton>
    <with-projection>false</with-projection>
    <!--<pression-explicite>true</pression-explicite>-->
    <schema-csts>0</schema-csts>
     <deltat-init>50.e-12</deltat-init>
     <deltat-min>1.e-14</deltat-min>
     <deltat-max>1.e-10</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>1.e-7</final-time>
    
    <time-history>
    <periode>3</periode>
    </time-history>
    
    <boundary-condition>
      <surface>XMIN</surface>
      <type>GeometricPressure</type>
      <value>10.e9</value>
      <Tdebut>1.e-9</Tdebut>   
      <Tfin>10.e-9</Tfin>      
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
