#! /usr/bin/awk -f
{
  sequence=$NF

  ls = length(sequence)
  is = 1
  fld  = 1
  while (fld < NF)
  {
     if (fld == 1){printf ">"}
     printf "%s " , $fld
     if (fld == NF-1){
        printf "\n"
      }
      fld = fld+1
  }
  while (is <= ls){
    printf "%s\n", substr(sequence,is,60)
    is=is+60
  }
}
