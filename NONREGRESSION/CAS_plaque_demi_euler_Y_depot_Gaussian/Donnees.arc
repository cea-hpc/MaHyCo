<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>CAS_lag_tache_focale</title>
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
      <variable>VelocityGradient</variable>
      <variable>DeformationRate</variable>
      <variable>StrainTensorXX</variable>
      <variable>PlasticDeformationVelocity</variable>
      <variable>PlasticDeformation</variable>
      <variable>SoundSpeed</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
    <output-history-period>1</output-history-period>
    <output-period>500</output-period>
  </arcane-post-processing>

  <mesh nb-ghostlayer="2" ghostlayer-builder-version="3">
    <meshgenerator>
    <cartesian>
       <nsd>1 4</nsd> 
       <origine>0.0 -.25e-3 0.0</origine>
       <lx nx='200' prx='1.0'>100.e-6</lx>
       <ly ny='1000' pry='1.0'>.5e-3</ly>
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
    <densite-initiale>2710.00</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <eos-model name="MieGruneisen">
      <adiabatic-cst>2.</adiabatic-cst>
      <rho0>2785.</rho0>
      <s1>1.339</s1>
      <s2>0.</s2>
      <s3>0.</s3>
      <gamma0>1.6</gamma0>
      <a>0.48</a>
      <c0>5328.</c0>
      <energie-ref>0.</energie-ref>
      <specific-heat>897.</specific-heat>
      <temperature-ref>300.</temperature-ref>
    </eos-model>
    <linear-pseudo-coeff>0.5</linear-pseudo-coeff>
    <quadratic-pseudo-coeff>1.2</quadratic-pseudo-coeff>
    <elasto-model name="NoModel">
    </elasto-model>
    <energy-depot>
    	<valeur-source-energie>3.7e16</valeur-source-energie>
	<type>DepotSuperGaussian</type>
 	<power>8.</power>
	<dependanceY>0.8e4</dependanceY>
     	<cutoffY>1.25e-4</cutoffY>
     	<cutoffX>10.e-6</cutoffX>
      	<Tdebut>0.</Tdebut>   
      	<Tfin>10.e-9</Tfin>
    </energy-depot>

  </environment>
  
   
   <cas-model name="ONECELL">
   <cas-test>0</cas-test>
   </cas-model>
    <remap name="RemapADI">
        <ordre-projection>2</ordre-projection>
    	<is-euler-scheme>Y</is-euler-scheme>
    </remap>
    <with-projection>true</with-projection>
    <pseudo-centree>0</pseudo-centree>
    <with-newton>false</with-newton>
    <pression-explicite>true</pression-explicite>
    <schema-csts>0</schema-csts>
    <deltat-init>1.e-13</deltat-init>
    <deltat-min>1.e-15</deltat-min>
    <deltat-max>1.e-13</deltat-max>
    <deltat-constant>false</deltat-constant>
    <longueur-caracteristique>monodimX</longueur-caracteristique>
     
    <final-time>1.e-11</final-time>
    
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
