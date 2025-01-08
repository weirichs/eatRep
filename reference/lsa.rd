<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no"><title>Achievement data from two large-scale assessments of 2010 and 2015. — lsa • eatRep</title><script src="../deps/jquery-3.6.0/jquery-3.6.0.min.js"></script><meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no"><link href="../deps/bootstrap-5.3.1/bootstrap.min.css" rel="stylesheet"><script src="../deps/bootstrap-5.3.1/bootstrap.bundle.min.js"></script><link href="../deps/font-awesome-6.5.2/css/all.min.css" rel="stylesheet"><link href="../deps/font-awesome-6.5.2/css/v4-shims.min.css" rel="stylesheet"><script src="../deps/headroom-0.11.0/headroom.min.js"></script><script src="../deps/headroom-0.11.0/jQuery.headroom.min.js"></script><script src="../deps/bootstrap-toc-1.0.1/bootstrap-toc.min.js"></script><script src="../deps/clipboard.js-2.0.11/clipboard.min.js"></script><script src="../deps/search-1.0.0/autocomplete.jquery.min.js"></script><script src="../deps/search-1.0.0/fuse.min.js"></script><script src="../deps/search-1.0.0/mark.min.js"></script><!-- pkgdown --><script src="../pkgdown.js"></script><meta property="og:title" content="Achievement data from two large-scale assessments of 2010 and 2015. — lsa"><meta name="description" content="This example data set contains fictional achievement scores of 11637 students from three countries
and two times of measurement in two domains (reading and listening comprehension) in the long format.
The data set contains nested multiple imputed plausible values of achievement scores as well as some
demographic variables. Illustrating trend analyses, data from two fictional time points (2010 and 2015)
are included.
The data set can be used for several illustration purposes. For example, if only multiple imputation
should be considered (without nesting), simply use only cases from the first nest (by subsetting). If
only one time of measurement should be considered (i.e., without any trend analyses), simply choose
only cases from 2010 or 2015. If only reading or listening should be considered, choose the desired
domain by subsetting according to the domain column."><meta property="og:description" content="This example data set contains fictional achievement scores of 11637 students from three countries
and two times of measurement in two domains (reading and listening comprehension) in the long format.
The data set contains nested multiple imputed plausible values of achievement scores as well as some
demographic variables. Illustrating trend analyses, data from two fictional time points (2010 and 2015)
are included.
The data set can be used for several illustration purposes. For example, if only multiple imputation
should be considered (without nesting), simply use only cases from the first nest (by subsetting). If
only one time of measurement should be considered (i.e., without any trend analyses), simply choose
only cases from 2010 or 2015. If only reading or listening should be considered, choose the desired
domain by subsetting according to the domain column."></head><body>
    <a href="#main" class="visually-hidden-focusable">Skip to contents</a>


    <nav class="navbar navbar-expand-lg fixed-top bg-light" data-bs-theme="light" aria-label="Site navigation"><div class="container">

    <a class="navbar-brand me-2" href="../index.html">eatRep</a>

    <small class="nav-text text-muted me-auto" data-bs-toggle="tooltip" data-bs-placement="bottom" title="">0.15.0</small>


    <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbar" aria-controls="navbar" aria-expanded="false" aria-label="Toggle navigation">
      <span class="navbar-toggler-icon"></span>
    </button>

    <div id="navbar" class="collapse navbar-collapse ms-3">
      <ul class="navbar-nav me-auto"><li class="nav-item"><a class="nav-link" href="../articles/eatRep.html">Get started</a></li>
<li class="active nav-item"><a class="nav-link" href="../reference/index.html">Reference</a></li>
<li class="nav-item"><a class="nav-link" href="../news/index.html">Changelog</a></li>
      </ul><ul class="navbar-nav"><li class="nav-item"><form class="form-inline" role="search">
 <input class="form-control" type="search" name="search-input" id="search-input" autocomplete="off" aria-label="Search site" placeholder="Search for" data-search-index="../search.json"></form></li>
<li class="nav-item"><a class="external-link nav-link" href="https://github.com/weirichs/eatRep/" aria-label="GitHub"><span class="fa fab fa-github fa-lg"></span></a></li>
      </ul></div>


  </div>
</nav><div class="container template-reference-topic">
<div class="row">
  <main id="main" class="col-md-9"><div class="page-header">

      <h1>Achievement data from two large-scale assessments of 2010 and 2015.</h1>

      <div class="d-none name"><code>lsa.rd</code></div>
    </div>

    <div class="ref-description section level2">
    <p>This example data set contains fictional achievement scores of 11637 students from three countries
and two times of measurement in two domains (reading and listening comprehension) in the long format.
The data set contains nested multiple imputed plausible values of achievement scores as well as some
demographic variables. Illustrating trend analyses, data from two fictional time points (2010 and 2015)
are included.</p>
<p>The data set can be used for several illustration purposes. For example, if only multiple imputation
should be considered (without nesting), simply use only cases from the first nest (by subsetting). If
only one time of measurement should be considered (i.e., without any trend analyses), simply choose
only cases from 2010 or 2015. If only reading or listening should be considered, choose the desired
domain by subsetting according to the <code>domain</code> column.</p>
    </div>

    <div class="section level2">
    <h2 id="ref-usage">Usage<a class="anchor" aria-label="anchor" href="#ref-usage"></a></h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="fu"><a href="https://rdrr.io/r/utils/data.html" class="external-link">data</a></span><span class="op">(</span><span class="va">lsa</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div class="section level2">
    <h2 id="format">Format<a class="anchor" aria-label="anchor" href="#format"></a></h2>
    <p>'data.frame':   77322 obs. of  25 variables</p><dl><dt>year</dt>
<dd><p>Year of evaluation</p></dd>

    <dt>idstud</dt>
<dd><p>individual student identification</p></dd>

    <dt>idclass</dt>
<dd><p>class identifier</p></dd>

    <dt>wgt</dt>
<dd><p>Total case weight</p></dd>

    <dt>L2wgt</dt>
<dd><p>School weight (level 2 weight)</p></dd>

    <dt>L1wgt</dt>
<dd><p>Student weight (level 1 weight)</p></dd>

    <dt>jkzone</dt>
<dd><p>jackknifing zone (jk2)</p></dd>

    <dt>jkrep</dt>
<dd><p>jackknife replicate</p></dd>

    <dt>imp</dt>
<dd><p>Number of imputation</p></dd>

    <dt>nest</dt>
<dd><p>Number of nest (for nested imputation only)</p></dd>

    <dt>country</dt>
<dd><p>The country an examinee stems from</p></dd>

    <dt>sex</dt>
<dd><p>student's sex</p></dd>

    <dt>ses</dt>
<dd><p>student's socio-economical status</p></dd>

    <dt>mig</dt>
<dd><p>student's migration background</p></dd>

  	<dt>domain</dt>
<dd><p>The domain the corresponding score belongs to</p></dd>

  	<dt>score</dt>
<dd><p>student's achievement score (corresponding to the domain reading or listening, and to the imputation 1, 2, or 3)</p></dd>

  	<dt>comp</dt>
<dd><p>student's competence level</p></dd>

  	<dt>failMin</dt>
<dd><p>dichotomous indicator whether the student fails to fulfill the minimal standard</p></dd>

  	<dt>passReg</dt>
<dd><p>dichotomous indicator whether the student fulfills at least the regular standard</p></dd>

  	<dt>passOpt</dt>
<dd><p>dichotomous indicator whether the student fulfills the optimal standard</p></dd>

  	<dt>leSore</dt>
<dd><p>linking error of each student's achievement score</p></dd>

  	<dt>leComp</dt>
<dd><p>linking error of each student's competence level</p></dd>

  	<dt>leFailMin</dt>
<dd><p>linking error of each student's indicator of failing to fulfill the minimal standard</p></dd>

  	<dt>lePassReg</dt>
<dd><p>linking error of each student's indicator of fulfilling the regular standard</p></dd>

  	<dt>lePassOpt</dt>
<dd><p>linking error of each student's indicator of fulfilling the optimal standard</p></dd>


</dl></div>
    <div class="section level2">
    <h2 id="source">Source<a class="anchor" aria-label="anchor" href="#source"></a></h2>
    <p>Simulated data</p>
    </div>

  </main><aside class="col-md-3"><nav id="toc" aria-label="Table of contents"><h2>On this page</h2>
    </nav></aside></div>


    <footer><div class="pkgdown-footer-left">
  <p>Developed by Sebastian Weirich, Martin Hecht, Karoline Sachse, Benjamin Becker.</p>
</div>

<div class="pkgdown-footer-right">
  <p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.1.1.</p>
</div>

    </footer></div>





  </body></html>

