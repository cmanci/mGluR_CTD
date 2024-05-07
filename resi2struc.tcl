#code to measure the secondary structure formed per residue in an IDP throughout an MD trajectory
#Dagan Marx, PhD
#output is a text file with 3 columns, for each type of secondary structure

set outfile [open avg_ss.dat w]
for {set j 261} {$j < 312} {incr j} {

  #set outfile [open avg_ss.dat w]
  set lookup {H G I}
  set lookup2 {T C}
  set lookup3 {E B}
  set frame_num [molinfo top get numframes]
  set full [atomselect top "protein and name CA and resid $j"]
  set len [llength [$full get resid]]
  $full delete
  set counterc 0
  set counterh 0
  set counterb 0
  set percenth 0
  set percentc 0
  set percentb 0
  for {set i 0} {$i < $frame_num} {incr i} {
      animate goto $i
      set sel [atomselect top "protein and name CA and resid $j"]
      mol ssrecalc top
      set struc_string [$sel get structure]
      #puts $struc_string
      set helix 0
      set beta 0
      set coil 0
      foreach letter $lookup {
          set temp1 [expr {[llength [split $struc_string $letter]] - 1}]
          #puts $temp
          incr helix $temp1
      }
      foreach letter $lookup2 {
          set temp2 [expr {[llength [split $struc_string $letter]] - 1}]
          #puts $temp
          incr coil $temp2
      }
      foreach letter $lookup3 {
          set temp3 [expr {[llength [split $struc_string $letter]] - 1}]
          #puts $temp
          incr beta $temp3
      }
      set percenth [expr {double($helix) / double($len) * 100}]
      set percentc [expr {double($coil) / double($len) * 100}]
      set percentb [expr {double($beta) / double($len) * 100}]
      #puts $outfile "$i\t$percenth\t$percentc\t$percentb"
      #puts $percentc
      set counterc [expr {$counterc + $percentc }]
    set counterh [expr {$counterh + $percenth }]
    set counterb [expr {$counterb + $percentb }]
    $sel delete
  }

  set avgc [expr $counterc/$frame_num]
  set avgh [expr $counterh/$frame_num]
  set avgb [expr $counterb/$frame_num]
  puts $outfile "$j\t$avgh\t$avgc\t$avgb"
}
close $outfile

#puts $avgc
