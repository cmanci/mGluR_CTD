#Code to measure position of a specific residue, or set of residues, relative to the phosphate plane of the membrane
#Dagan Marx, PhD
#produces a text file containing the position relative to the phosphate time as a function of timestep in your VMD trajectory
#important to import .dcd files aligned to the membrane and not the protein so the relative positions of the phosphate plane are set such that one leaflet has a positive z-coordinate value and one leaflet has a negative z-coordinate value

#choose residue number or set of residues using VMD selection algebra
set residue "resid 311"
set phos_upper [atomselect top "name P and  z> 0"]
set phos_lower [atomselect top "name P and z< 0"]

#select which parts of the residue you want, i.e. sidechain or whole amino acid
set residue_ca [atomselect top "protein and $residue"]

#set name of outfile 
set outfile2 [open CtermCTD_posrelphos.txt w]
set nf [molinfo top get numframes]


for {set i 0} {$i < $nf} {incr i} {
  $phos_upper frame $i
  $phos_lower frame $i
  $residue_ca frame $i
  set av_phosz_upper [lindex [measure center $phos_upper weight mass] 2]
  set av_phosz_lower [lindex [measure center $phos_lower weight mass] 2]
  set ca_z [lindex [measure center $residue_ca weight mass] 2]
  set phos_center [expr (($av_phosz_upper + $av_phosz_lower)/2)] 
  if {$ca_z > $phos_center} { set ca_adjz [expr $ca_z - $av_phosz_upper] }
  if {$ca_z < $phos_center} { set ca_adjz [expr $ca_z - $av_phosz_lower] }
  set ca_adjz [expr $ca_z - $av_phosz_lower] 
  set ca_adjz [expr $ca_adjz * -1]
  puts $outfile2 $ca_adjz

}


close $outfile2
