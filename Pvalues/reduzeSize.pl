while(<>){
$i++;
if ($i<140 or $i>143){
print
}else{
@t=split(/ L /);
print $t[0];
%h=();
for $i(1..$#t){
@t1=split(/ /,$t[$i]);
next if index($t1[0],'.')<0 or index($t1[1],'.')<0;
$x=substr($t1[0],0,index($t1[0],'.')+2);
$y=substr($t1[1],0,index($t1[1],'.')+2);
print " L $x $y\n" unless defined $h{"$x $y"};
$h{"$x $y"}=1
}
print " L $t[$#t]";
}
}
