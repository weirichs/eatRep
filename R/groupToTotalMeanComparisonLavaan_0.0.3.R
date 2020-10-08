# groupToTotalMeanComparisonLavaan, 0.0.1 (2020-10-05)

# Notes:
#  -- complete data set is assumed (no missings)
#  -- lavaan seems to go by order of groups in group variable in data set, not by defined levels of factor!
#  -- lavaan seems to ignore sampling.weights=weight.var for Poisson model, only means (from regression) are calculated with sample weights
#     => total mean calculation needs to be adjusted (for stochastic group sizes)

# Arguments:
#  d                     ... data frame (no missing values)
#  y.var                 ... name of y variable
#  group.var             ... name of group variable
#  weight.var            ... name of sample weights variable
#  heterogeneous         ... heterogeneous group variances
#  stchgrsz              ... stochastic group sizes
#  extended.results      ... if FALSE only group to total mean differences are returned; # if TRUE Mtotal is additionally included in returned results
#  lavaan.syntax.path    ... path to lavaan syntax (just for reference purposes, not needed)
#  run                   ... if TRUE lavaan is run and results are returned, if FALSE lavaan model syntax is returned
#  lavaan.summary.output ... show lavaan's summary output on console

# Value:
#  depending on 'run' either a data frame with results (ordered by factor levels or by occurence in data frame) or lavaan syntax as character string

groupToTotalMeanComparisonLavaan <- function( d, y.var, group.var, weight.var=NULL, heterogeneous=FALSE, stchgrsz=FALSE, extended.results=FALSE, lavaan.syntax.path=NULL, run=TRUE, lavaan.summary.output=FALSE ){

		# lavaan, janitor is required
#		require( lavaan )
#		require( janitor )

		# lavaanVersion
#		lavaanVersion <- installed.packages()["lavaan","Version"]
		# first, second, third digit of lavaanVersion
#		fd <- suppressWarnings( as.numeric( sub( "^(\\d+)\\..*$", "\\1", lavaanVersion ) ) )
#		sd <- suppressWarnings( as.numeric( sub( "^\\d+\\.(\\d+).*$", "\\1", lavaanVersion ) ) )
#		td <- suppressWarnings( as.numeric( sub( "^\\d+\\.\\d+-(\\d).*$", "\\1", lavaanVersion ) ) )
#		if( !any( is.na( c( fd, sd, td ) ) ) ){
			# version as numeric
#			v <- as.numeric( paste0( "1", formatC( fd, width=3, flag="0" ), "1", formatC( sd, width=3, flag="0" ), "1", formatC( td, width=3, flag="0" ) ) )
			# stop if not at least lavaan version 0.6-7
#			if( v < 100010061007 ) stop( "lavaan version 0.6-7 or newer is required" )
#		}

		# variable names inputted by user
		vars <- c( y.var, group.var, weight.var )

		# check if in colnames of data frame
		if( !all( indfr <- vars %in% colnames( d ) ) ) stop( paste0( "variable(s) ", paste(paste0("'",vars[!indfr],"'"),collapse=", "), " not in data" ) )

		# check for missings
		if( any( sapply( vars, function(var) any( is.na( d[, var] ) ) ) ) ) stop( "data contains missings, not supported" )

		# if factor, keep levels (but only existing ones in the data frame)
		if( is.factor( d[,group.var] ) ) { orig.levels <- levels( d[,group.var] ); orig.levels <- orig.levels[orig.levels %in% unique( as.character( d[,group.var] ))] } else { orig.levels <- NULL }

		# convert group variable to character
		d[,group.var] <- as.character ( d[,group.var] )

		# original group labels
		orig.groups <- unique( d[,group.var] )

		# clean group labels
		group.lab <- make_clean_names( orig.groups )

		# forbidden group labels
		if( any( group.lab %in% (forbidden.gr.lab <- c("Mtotal", "N")) ) ) stop( paste0( "group labels must not be ", paste( forbidden.gr.lab, collapse=" | " ) ) )

		# relative group frequencies (fixed, i.e., from sample), relevant for stchgrsz=FALSE
		if( is.null( weight.var ) ){
				rgf <- sapply( orig.groups, function( g ) length( d[d[,group.var] %in% g,y.var] ) )
				# just for reference purpose, calculate group means and total mean (it's outputted into the lavaan syntax)
				Mg <- sapply( orig.groups, function( g ) mean( d[d[,group.var] %in% g,y.var] ) )
				Mtotal <- mean( d[,y.var] )
		} else {
				rgf <- sapply( orig.groups, function( g ) sum( d[d[,group.var] %in% g,weight.var] ) )
				# just for reference purpose, calculate weighted group means and weighted total mean (it's outputted into the lavaan syntax)
				Mg <- sapply( orig.groups, function( g ) weighted.mean( d[d[,group.var] %in% g,y.var], d[d[,group.var] %in% g,weight.var] ) )
				Mtotal <- weighted.mean( d[,y.var], d[,weight.var] )
		}
		names( rgf ) <- group.lab
		N <- sum( rgf )
		relfreq <- rgf / N
		names( Mg ) <- group.lab

		# "weighted by"-string (just for comments in lavaan code)
		w.by1 <- ifelse(!is.null(weight.var),paste0( "\n# weighted by sample weights variable '", weight.var, "' via laavan argument sampling.weights='",weight.var,"'\n" ), "\n" )
		w.by2 <- ifelse(!is.null(weight.var) & !stchgrsz,paste0( "\n# weighted by using the sample weights variable '", weight.var, "'" ), "" )

		# lavaan model
		m <- paste0( 
				'# dependent variable is regressed on intercept in each group (gives group means)',w.by1,
				y.var,' ~ c(',paste(paste0("M",group.lab),collapse=","),')*1\n',
				ifelse( !heterogeneous, paste0( '\n# homogeneous group variance\n',y.var,' ~~ c(', paste(rep("Var",length(group.lab)),collapse=","), ")*",y.var, "\n\n" ), "\n" ),
				ifelse( stchgrsz, paste0('# relative group frequencies (Poisson model)',w.by1,'group % c(',paste(paste0('W', group.lab ),collapse=','),')*w\nN := ',paste(paste0('exp( W',group.lab ," )"),collapse=" + "),'\n',paste(paste0("relfreq", group.lab ," := exp( W",group.lab," )/N" ),collapse="\n"), "\n" ), paste0( "# fixed group sizes (i.e., from sample)", w.by2, "\n", paste( paste0( "relfreq", names(relfreq), " := ", relfreq ), collapse="\n" ), "\n" ) ),
				'\n# total mean',
				'\nMtotal := ',paste0( paste( paste0( "relfreq", group.lab, "*M", group.lab ), collapse="+" ), "\n" ),
				'\n# group mean to total mean difference (sample values behind #)\n',
				paste(paste0(group.lab," := M",group.lab," - Mtotal # ",Mg, " - ", Mtotal, " = ", Mg-Mtotal ),collapse="\n")
		)

		# write lavaan syntax to path
		if( !is.null( lavaan.syntax.path ) ){
				tried <- try( write( m, lavaan.syntax.path ) )
				if( inherits( tried, "try-error" ) ) warning( "Couldn't create lavaan syntax" )
		}
			
		if( run ){

				# run lavaan model
				r <- sem( m, data=d, group=group.var, sampling.weights=weight.var, estimator="MLR" )

				# results
				if( !lavaan.summary.output ) sink("nul")
				smr <- summary( r )
				if( !lavaan.summary.output ) sink()

				# results data frame
				res <- smr[[1]]
				colnames( res )[colnames( res ) %in% "label"] <- "par"

				# keep only relevant results
				res2 <- res[ res$par %in% group.lab, rel.cols <- c("par","est","se","z","pvalue")]

				# restore parameter labels (as before cleaning)
				res3 <- res2
				res3$par <- orig.groups

				# if group variable was originally factor, sort results by levels
				if ( !is.null( orig.levels ) ){
					res3 <- res3[ match( orig.levels, res3$par ), ]
				}

				# add extended results (Mtotal)
				if( extended.results ) res3 <- rbind( res3, res[ res$par %in% "Mtotal", rel.cols ] )

				# rownames
				rownames( res3 ) <- seq( along=rownames( res3 ) )

				# return
				return( res3 )
				
		} else { return( m ) }
}
