<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_Lag_ContactHerz3D</title>
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
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
    <output-history-period>500</output-history-period>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>2 2 1</nsd> 
       <origine>-50.e-6 -50.e-6 0.</origine>
       <lx nx='50' prx='1.0'>100.e-6</lx>
       <ly ny='50' pry='1.0'>100.e-6</ly>
       <lz nz='10' prz='1.0'>50.e-6</lz>
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
    <borne-sup>91.e-6 6.e-6 0.</borne-sup>
    <borne-inf>90.e-6 5.e-6 0.</borne-inf>
    </bilan>
  </time-history>

  <!-- Configuration du module hydrodynamique -->
  <mahyco>
   <material><name>CU_mat</name></material>
  <environment>
    <name>CU</name>
    <material>CU_mat</material>
    <densite-initiale>8930.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <eos-model name="MieGruneisen">
      <adiabatic-cst>2.02</adiabatic-cst>
      <rho0>8930.</rho0>
      <s1>1.489</s1>
      <s2>0.</s2>
      <s3>0.</s3>
      <gamma0>2.02</gamma0>
      <a>0.47</a>
      <c0>3940.</c0>
      <energie-ref>0.</energie-ref>
      <specific-heat>385.</specific-heat>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
    <elasto-model name="DefaultModel">
      <yandg-model name="SCG">
      <Y0>1.2e8</Y0>
      <Beta>36</Beta>
      <n>0.45</n>
      <Epsilon_init>0.</Epsilon_init>
      <Mu0>4.77e10</Mu0>
      <GPP>1.3356</GPP>
      <GPT>-1.8126e+7</GPT>
      <Ymax>6.4e8</Ymax>
      <Eint_fus>460e3</Eint_fus>
      </yandg-model>
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
     
    <final-time>2.e-9</final-time>
    
    <boundary-condition>
      <surface>ZMAX</surface>
      <type>ContactHerzPressure</type>
      <value>5.e9</value>
      <dependanceX>100.e-12</dependanceX> 
      <dependanceY>200.e-12</dependanceY> 
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
