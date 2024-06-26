<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_lag_Chauffage_isochore</title>
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
    <output-history-period>100</output-history-period>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>1 1</nsd> 
       <origine>-5.e-6 -5.e-6</origine>
       <lx nx='5' prx='1.0'>10.e-6</lx>
       <ly ny='5' pry='1.0'>10.e-6</ly>
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
  <material><name>SN_mat</name></material>
  <environment>
    <name>SN</name>
    <material>SN_mat</material>
    <densite-initiale>8200.</densite-initiale>
    <energie-initiale>7.738376</energie-initiale>
    <eos-model name="Autre">
      <fichier-coeff>ee.CineTest22#.Sn.00#.coeff</fichier-coeff>
    </eos-model>
    <energy-depot>
    <type>DepotConstant</type>
    <valeur-source-energie>13.2e14</valeur-source-energie>
    </energy-depot>
  </environment>
  
   
   <cas-model name="ONECELL">
   <cas-test>0</cas-test>
   </cas-model>
   
    <pseudo-centree>0</pseudo-centree>
    <with-newton>false</with-newton>
    <with-projection>false</with-projection>
    <pression-explicite>true</pression-explicite>
    <schema-csts>0</schema-csts>
     <deltat-init>1.e-13</deltat-init>
     <deltat-min>1.e-14</deltat-min>
     <deltat-max>1.e-13</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>5.e-10</final-time>

   <time-history>
    <periode>10000</periode>
    <borne-inf>0.e-6 0.e-6 -1.</borne-inf>
    <borne-sup>4.e-6 4.e-6 1.</borne-sup>
    </time-history>

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
		
  </mahyco>
</case>
