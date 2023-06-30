<?xml version='1.0'?>
<case codeversion="1.0" codename="Mahyco" xml:lang="en">
  <arcane>
    <title>Exemple Arcane d'un module Hydro tres simplifie</title>
    <timeloop>MahycoLoop</timeloop>
  </arcane>

  <arcane-post-processing>
    <output-period>50</output-period>
    <output>
      <variable>CellMass</variable>
      <variable>Pressure</variable>
      <variable>Density</variable>
      <variable>Velocity</variable>
      <variable>NodeMass</variable>
      <variable>InternalEnergy</variable>
      <variable>PseudoViscosity</variable>
      <variable>Materiau</variable>
      <variable>FracPhase1</variable>
      <variable>FracPhase2</variable>
      <variable>FracPhase3</variable>
    </output>
    <format>
      <binary-file>false</binary-file>
    </format>
    <output-history-period>100</output-history-period>
  </arcane-post-processing>

  <mesh>
    <meshgenerator>
    <cartesian>
       <nsd>1 1 1</nsd> 
       <origine>0.0 0.0 0.0</origine>
       <lx nx='4' prx='1.0'>10.e-6</lx>

       <ly ny='4' pry='1.0'>10.e-6</ly>

       <lz nz='1' prz='1.0'>1.e-6</lz>
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
  <material><name>SN_mat</name></material>
  <environment>
    <name>SN</name>
    <material>SN_mat</material>
    <densite-initiale>7286.62</densite-initiale>
    <pression-initiale>1.e5</pression-initiale>
    <energie-initiale>7.738376</energie-initiale>
    <eos-model name="Autre">
    <adiabatic-cst>13.2e15</adiabatic-cst>
    </eos-model> 
  </environment>
  
   
   <cas-model name="ONECELL">
   <cas-test>0</cas-test>
   </cas-model>
   
    <pseudo-centree>0</pseudo-centree>
    <with-newton>false</with-newton>
    <energy-deposit>true</energy-deposit>
    <with-projection>false</with-projection>
    <pression-explicite>true</pression-explicite>
    <schema-csts>0</schema-csts>
     <deltat-init>1.e-13</deltat-init>
     <deltat-min>1.e-14</deltat-min>
     <deltat-max>1.e-13</deltat-max>
    <longueur-caracteristique>racine-cubique-volume</longueur-caracteristique>
     
    <final-time>4.7e-12</final-time>

    <time-history>
    <periode>3</periode>
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
    <boundary-condition>
      <surface>ZMIN</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
    <boundary-condition>
      <surface>ZMAX</surface>
      <type>Vz</type>
      <value>0.</value>
    </boundary-condition>
		
  </mahyco>
</case>
