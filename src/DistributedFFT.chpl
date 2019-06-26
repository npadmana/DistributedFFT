/* Documentation for DistributedFFT */
prototype module DistributedFFT {

  use BlockDist;
  use FFTW;
  use FFTW.C_FFTW;


  /* Convenience constructor for grids */
  proc newSlabDom(dom: domain) where isRectangularDom(dom) {
    if dom.rank < 2 then compilerError("The domain must be at least 2D");
    if ((dom.dim(1).size)%numLocales !=0) then halt("numLocales must divide the first dimension");
    const targetLocales = reshape(Locales, {0.. #numLocales, 0..0, 0..0});
    return dom dmapped Block(boundingBox=dom, targetLocales=targetLocales);
  }

  proc newSlabDom(sz) where isHomogeneousTupleType(sz.type) {
    var tup : (sz.size)*range;
    for param ii in 1..sz.size do tup(ii) = 0.. #sz(ii);
    return newSlabDom({(...tup)});
  }


}