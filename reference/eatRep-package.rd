<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Statistical analyses in complex survey designs with multiple imputed data and trend estimation. — eatRep-package • eatRep</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js" integrity="sha512-v2CJ7UaYy4JwqLDIrZUI/4hqeoQieOmAZNXBeQyjo21dadnwR+8ZaIJVT8EE2iyI61OV8e6M8PP2/4hpQINQ/g==" crossorigin="anonymous" referrerpolicy="no-referrer"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Statistical analyses in complex survey designs with multiple imputed data and trend estimation. — eatRep-package"><meta property="og:description" content="The package provide functions to computes some basic statistic operations—(adjusted) means, standard deviations,
  frequency tables, percentiles and generalized linear models—in complex survey designs comprising multiple
  imputed variables and/or a clustered sampling structure which both deserve special procedures at least in
  estimating standard errors. In large-scale assessments, standard errors are comprised of three components:
  the measurement error, the sampling error, and (if trend estimation of at least two times of measurement
  are involved) the linking error.
Measurement error: In complex surveys or large-scale assessments, measurement errors are taken
  into account by the mean of multiple imputed variables. The computation of standard errors for the mean
  of a multiple imputed variable (e.g. plausible values) involves the formulas provided by Rubin (1987).
  Computing standard errors for the mean of a nested imputed variable involves the formulas provided by
  Rubin (2003). Both methods are implemented in the package. The estimation of \(R^2\) and adjusted
  \(R^2\) in linear and generalized linear regression models with multiple imputed data sets is
  realized using the methods provided in Harel (2009).
Sampling error: Computation of sampling errors of variables which stem from a clustered design may
  involve replication methods like balanced repeated replicate (BRR), bootstrap or Jackknife methods.
  See Westat (2000), Foy, Galia &amp;amp; Li (2008), Rust and Rao (1996), and Wolter (1985) for details. To date,
  the Jackknife-1 (JK1), Jackknife-2 (JK2) and the Balanced Repeated Replicates (BRR; optionally with Fay's
  method) procedures are supported.
Linking error: Lastly, standard errors for trend estimates may involve incorporating
  linking errors to account for potential differential item functioning or item parameter drift.
  eatRep allows to account for linking error when computing standard errors for trend
  estimates. Standard error estimation is conducted according to the operational practice in
  PISA, see equation 5 in Sachse &amp;amp; Haag (2017).
The package eatRep is designed to combine one or several error types which is necessary,
  for example, if (nested) multiple imputed data are used in clustered designs. Considering the
  structure is relevant especially for the estimation of standard errors. The estimation of national
  trends requires a sequential analysis for both measurements and a comparison of estimates between them.
Technically, eatRep is a wrapper for the survey package (Lumley, 2004). Each function in
  eatRep corresponds to a specific function in survey which is called repeatedly during the analysis.
  Hence, a nested loop is used. We use &amp;#8220;trend replicates&amp;#8221; in the outer loop, &amp;#8220;imputation replicates&amp;#8221;
  in the middle loop to account for multiple imputed data, and &amp;#8220;cluster replicates&amp;#8221; in the inner loop to
  account for the clustered sampling structure. While the functional principle of survey is based on
  replication of standard analyses, eatRep is based on replication of survey analyses to take
  multiple imputed data into account. More recent versions of the package additionally allow estimations using
  the BIFIEsurvey package instead of survey which provide substantial advantages in terms of speed.
For each imputed data set in each measurement, i.e. in the inner loop, the eatRep function first creates
  replicate weights based on the primary sampling unit (PSU) variable and the replication indicator variable. In
  the jackknife procedure, the first one is often referred to as &amp;#8220;jackknife zone&amp;#8221;, whereas the second one
  is often referred to as &amp;#8220;jackknife replicate&amp;#8221;. The number of distinct units in the PSU variable defines
  the number of replications which are necessary due to the clustered structure. A design object is created and
  the appropriate survey function is called. The process is repeated for each imputed dataset and the
  results of the analyses are pooled. The pooling procedure varies in relation to the type of variable to be
  pooled. For examples, means or regression coefficients are pooled according to Rubin (1987) or Rubin (2003).
  \(R^2\) is pooled according to Harel (2009), using a Fisher z-transformation. Chi-square distributed values
  are pooled according to Thomas and Rao (1990) for clustered data and according to Enders (2010) and
  Allison (2002) for multiple imputed data. For trend analyses, the whole process is repeated two times
  (according to the two measurements) and the difference of the estimates are computed along with their
  pooled standard errors.
Without trend estimation, the outer loop has only one cycle (instead of two). Without multiple imputations,
  the middle loop has only one cycle. Without a clustered sampling structure (i.e, in a random sample), the
  inner loop has only one cycle. Without trend, imputation and clustered structure, no replication is performed
  at all. To compute simple mean estimates, for example, eatRep then simply calls mean instead
  of svymean from the survey package. A special case occurs with nested multiple imputation.
  We then have four loops in a nested structure. Hence, the corresponding analyses may take considerably
  computational effort.
Important note: Starting with version 0.10.0, several methods for the standard error estimation
  of cross level differences are implemented. Prior to version 0.10.0, the standard error for the difference
  between one single group (e.g., Belgium) and the total population (which is comprised of several states including
  Belgium) was estimated as if both groups would have been independent from each other. The standard errors,
  however, are biased then. Two new methods are now applicable using the argument crossDiffSE in
  repMean and provide unbiased standard errors—weighted effect coding (wec) and replication
  methods (rep); see, for example te Grotenhuis et al. (2017) and Weirich et al. (2021). The old method is still available by
  using crossDiffSE = &quot;old&quot;. Note that the default method now is weighted effect coding.
Second important note: Starting with version 0.13.0, function names have been changed due to
  inconsistent former denomination: Function jk2.mean now goes under the name of repMean,
  jk2.table was  renamed to repTable, jk2.quantile was  renamed to repQuantile,
  and jk2.glm now goes under the name of repGlm. The old functions are deprecated and will
  be removed in further package publications. Renaming was driven by the fact that the corresponding
  functions now have broader range of methods than only jackknife-2.
Third important note: Starting with version 0.15.0, the reporting function report was deprecated
  due to inefficient and error-prone programming. The new reporting function report2 has a new output
  format which provides an interface for the eatPlot package. Old functionality was roughly supplied using
  report, but if the 1:1 output of former version is requested, please use version 0.14.7."><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">


    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">eatRep</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">0.15.1</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li>
  <a href="../articles/eatRep.html">Get started</a>
</li>
<li>
  <a href="../reference/index.html">Reference</a>
</li>
<li>
  <a href="../news/index.html">Changelog</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/weirichs/eatRep/" class="external-link">
    <span class="fab fa-github fa-lg"></span>

  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->



      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Statistical analyses in complex survey designs with multiple imputed data and trend estimation.</h1>

    <div class="hidden name"><code>eatRep-package.rd</code></div>
    </div>

    <div class="ref-description">
    <p>The package provide functions to computes some basic statistic operations—(adjusted) means, standard deviations,
  frequency tables, percentiles and generalized linear models—in complex survey designs comprising multiple
  imputed variables and/or a clustered sampling structure which both deserve special procedures at least in
  estimating standard errors. In large-scale assessments, standard errors are comprised of three components:
  the measurement error, the sampling error, and (if trend estimation of at least two times of measurement
  are involved) the linking error.</p>
<p><strong>Measurement error:</strong> In complex surveys or large-scale assessments, measurement errors are taken
  into account by the mean of multiple imputed variables. The computation of standard errors for the mean
  of a multiple imputed variable (e.g. plausible values) involves the formulas provided by Rubin (1987).
  Computing standard errors for the mean of a nested imputed variable involves the formulas provided by
  Rubin (2003). Both methods are implemented in the package. The estimation of \(R^2\) and adjusted
  \(R^2\) in linear and generalized linear regression models with multiple imputed data sets is
  realized using the methods provided in Harel (2009).</p>
<p><strong>Sampling error:</strong> Computation of sampling errors of variables which stem from a clustered design may
  involve replication methods like balanced repeated replicate (BRR), bootstrap or Jackknife methods.
  See Westat (2000), Foy, Galia &amp; Li (2008), Rust and Rao (1996), and Wolter (1985) for details. To date,
  the Jackknife-1 (JK1), Jackknife-2 (JK2) and the Balanced Repeated Replicates (BRR; optionally with Fay's
  method) procedures are supported.</p>
<p><strong>Linking error:</strong> Lastly, standard errors for trend estimates may involve incorporating
  linking errors to account for potential differential item functioning or item parameter drift.
  <code>eatRep</code> allows to account for linking error when computing standard errors for trend
  estimates. Standard error estimation is conducted according to the operational practice in
  PISA, see equation 5 in Sachse &amp; Haag (2017).</p>
<p>The package <code>eatRep</code> is designed to combine one or several error types which is necessary,
  for example, if (nested) multiple imputed data are used in clustered designs. Considering the
  structure is relevant especially for the estimation of standard errors. The estimation of national
  trends requires a sequential analysis for both measurements and a comparison of estimates between them.</p>
<p>Technically, <code>eatRep</code> is a wrapper for the <code>survey</code> package (Lumley, 2004). Each function in
  <code>eatRep</code> corresponds to a specific function in <code>survey</code> which is called repeatedly during the analysis.
  Hence, a nested loop is used. We use “trend replicates” in the outer loop, “imputation replicates”
  in the middle loop to account for multiple imputed data, and “cluster replicates” in the inner loop to
  account for the clustered sampling structure. While the functional principle of <code>survey</code> is based on
  replication of standard analyses, <code>eatRep</code> is based on replication of <code>survey</code> analyses to take
  multiple imputed data into account. More recent versions of the package additionally allow estimations using
  the <code>BIFIEsurvey</code> package instead of <code>survey</code> which provide substantial advantages in terms of speed.</p>
<p>For each imputed data set in each measurement, i.e. in the inner loop, the <code>eatRep</code> function first creates
  replicate weights based on the primary sampling unit (PSU) variable and the replication indicator variable. In
  the jackknife procedure, the first one is often referred to as “jackknife zone”, whereas the second one
  is often referred to as “jackknife replicate”. The number of distinct units in the PSU variable defines
  the number of replications which are necessary due to the clustered structure. A design object is created and
  the appropriate <code>survey</code> function is called. The process is repeated for each imputed dataset and the
  results of the analyses are pooled. The pooling procedure varies in relation to the type of variable to be
  pooled. For examples, means or regression coefficients are pooled according to Rubin (1987) or Rubin (2003).
  \(R^2\) is pooled according to Harel (2009), using a Fisher <em>z</em>-transformation. Chi-square distributed values
  are pooled according to Thomas and Rao (1990) for clustered data and according to Enders (2010) and
  Allison (2002) for multiple imputed data. For trend analyses, the whole process is repeated two times
  (according to the two measurements) and the difference of the estimates are computed along with their
  pooled standard errors.</p>
<p>Without trend estimation, the outer loop has only one cycle (instead of two). Without multiple imputations,
  the middle loop has only one cycle. Without a clustered sampling structure (i.e, in a random sample), the
  inner loop has only one cycle. Without trend, imputation and clustered structure, no replication is performed
  at all. To compute simple mean estimates, for example, <code>eatRep</code> then simply calls <code>mean</code> instead
  of <code>svymean</code> from the <code>survey</code> package. A special case occurs with nested multiple imputation.
  We then have four loops in a nested structure. Hence, the corresponding analyses may take considerably
  computational effort.</p>
<p><em>Important note:</em> Starting with version 0.10.0, several methods for the standard error estimation
  of cross level differences are implemented. Prior to version 0.10.0, the standard error for the difference
  between one single group (e.g., Belgium) and the total population (which is comprised of several states including
  Belgium) was estimated as if both groups would have been independent from each other. The standard errors,
  however, are biased then. Two new methods are now applicable using the argument <code>crossDiffSE</code> in
  <code><a href="repMean.html">repMean</a></code> and provide unbiased standard errors—weighted effect coding (wec) and replication
  methods (rep); see, for example te Grotenhuis et al. (2017) and Weirich et al. (2021). The old method is still available by
  using <code>crossDiffSE = "old"</code>. Note that the default method now is weighted effect coding.</p>
<p><em>Second important note:</em> Starting with version 0.13.0, function names have been changed due to
  inconsistent former denomination: Function <code>jk2.mean</code> now goes under the name of <code><a href="repMean.html">repMean</a></code>,
  <code>jk2.table</code> was  renamed to <code><a href="repTable.html">repTable</a></code>, <code>jk2.quantile</code> was  renamed to <code><a href="repQuantile.html">repQuantile</a></code>,
  and <code>jk2.glm</code> now goes under the name of <code><a href="repGlm.html">repGlm</a></code>. The old functions are deprecated and will
  be removed in further package publications. Renaming was driven by the fact that the corresponding
  functions now have broader range of methods than only jackknife-2.</p>
<p><em>Third important note:</em> Starting with version 0.15.0, the reporting function <code><a href="report.html">report</a></code> was deprecated
  due to inefficient and error-prone programming. The new reporting function <code><a href="report.html">report2</a></code> has a new output
  format which provides an interface for the <code>eatPlot</code> package. Old functionality was roughly supplied using
  <code><a href="report.html">report</a></code>, but if the 1:1 output of former version is requested, please use version 0.14.7.</p>
    </div>


    <div id="details">
    <h2>Details</h2>

<table class="table table"><tr><td>Package:</td><td>eatRep</td></tr><tr><td>Type:</td><td>Package</td></tr><tr><td>Version:</td><td>0.15.1</td></tr><tr><td>Date:</td><td>2025-02-09</td></tr><tr><td>License:</td><td>GPL(&gt;=2)</td></tr></table></div>
    <div id="author">
    <h2>Author</h2>
    <p>Authors: Sebastian Weirich &lt;sebastian.weirich@iqb.hu-berlin.de&gt;, Martin Hecht &lt;martin.hecht@hu-berlin.de&gt;,
    Benjamin Becker &lt;b.becker@iqb.hu-berlin.de&gt;</p>
    </div>
    <div id="references">
    <h2>References</h2>
    <p>Allison, P. D. (2002). Missing data. Newbury Park, CA: Sage.</p>
<p>Enders, C. K. (2010). Applied missing data analysis. Guilford Press.</p>
<p>Foy, P., Galia , J. &amp; Li, I. (2008). Scaling the data from the TIMSS 2007 mathematics
  and science assessment. In J. F. Olson, M. O. Martin &amp; I. V. S. Mullis (ed.),
  <em>TIMSS 2007 Technical Report</em> (S. 225–280). Chestnut Hill, MA: TIMSS &amp; PIRLS
  International Study Center, Lynch School of Education, Boston College.</p>
<p>Harel, O. (2009): The estimation of \(R^2\) and adjusted \(R^2\) in incomplete data
  sets using multiple imputation. <em>Journal of Applied Statistics.</em> <b>36, 10</b>, 1109–1118.</p>
<p>Lumley, T. (2004). Analysis of complex survey samples. <em>Journal of Statistical Software</em> <b>9(1)</b>: 1–19</p>
<p>Rubin, D. B. (1987). <em>Multiple imputation for nonresponse in surveys.</em> New York: Wiley.</p>
<p>Rubin, D.B. (2003): Nested multiple imputation of NMES via partially incompatible MCMC.
  <em>Statistica Neerlandica</em> <b>57, 1</b>, 3–18.</p>
<p>Rust, K., &amp; Rao, JNK. (1996): Variance estimation for complex surveys using
  replication techniques. <em>Statistical Methods in Medical Research</em> <b>5</b>, 283–310.</p>
<p>Sachse, K. A. &amp; Haag, N. (2017). Standard errors for national trends in international
  large-scale assessments in the case of cross-national differential item functioning. <em>Applied
  Measurement in Education, 30</em>, (2), 102-116. http://dx.doi.org/10.1080/08957347.2017.1283315</p>
<p>Satorra, A., &amp; Bentler, P. M. (1994). Corrections to test statistics
		and standard errors in covariance structure analysis.</p>
<p>te Grotenhuis, M., Pelzer, B., Eisinga, R., Nieuwenhuis, R., Schmidt-Catran, A., &amp; Konig, R. (2017).
  When size matters: advantages of weighted effect coding in observational studies.
  <em>International Journal of Public Health.</em> <b>62</b>, 163–167.</p>
<p>Thomas, D. R. &amp; Rao, JNK (1990): Small-sample comparison of level and power for simple goodness-of-
  fit statistics under cluster sampling. JASA 82:630-636</p>
<p>Weirich, S., Hecht, M., Becker, B. et al. (2021). Comparing group means with the total mean in random samples,
  surveys, and large-scale assessments: A tutorial and software illustration. Behavior Research Methods.
  https://doi.org/10.3758/s13428-021-01553-1</p>
<p>Westat (2000). <em>WesVar.</em> Rockville, MD: Westat.</p>
<p>Wolter, K. M. (1985). <em>Introduction to variance estimation.</em> New York: Springer.</p>
    </div>

  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by Sebastian Weirich, Martin Hecht, Karoline Sachse, Benjamin Becker.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.1.1.</p>
</div>

      </footer></div>






  </body></html>

