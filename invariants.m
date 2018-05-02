//Compute the Molien series, the predicted degrees for the invariants (primary and secondary) and computes the primary, secondary, and fundamental invariants.

intrinsic Tesfct
 (nki::[RngIntElt],
  outputfile::Str,logfile::Str) ->
 RngIntElt
{Script for computing invariants}
  SetOutputFile(outputfile: Overwrite := true);
  SetLogFile(logfile: Overwrite := true); 
  //Attach("pack.m");
  K := RationalField();
  //G := MolSymGen(nki);
  //M<t> := MolienSeries(G);
  //M;
  M := 2;
  Write(outputfile, M);
return 0
end intrinsic;
