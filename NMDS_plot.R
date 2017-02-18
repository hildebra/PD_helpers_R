#plotting routines to plot vegan NMDS

plotNMDS_EF = function(mds,ef,class,Opt,IndFit=TRUE,maxTake=3,colObj=NULL,
		FeatNam=NULL,MetaInfo=NULL,usePoints = 0,groups=NULL,
		insertMetaToSit=FALSE,insertSitToSpecP=TRUE,ordiMethod="NMDS",
		clusObj=NULL,externalFacHull=NULL,justSites=F,FeatureDistrPlot=FALSE,
		TimeSeries=NULL,switchAxes=FALSE,hullType=NULL,ordiSpCol = NULL,
		displaySubs=NULL,doXtraLegend=NULL){
	#browser();
	#pvf=ef$factors$pvals;  names(pvf)=names(ef$factors$r); pvf=sort(pvf)
	totDim = mds$ndim; combis = t(combn(totDim,2)); Ncomb = dim(combis)[1]
	if (switchAxes){ combis = combis[,c(2,1),drop=F]}
	xlbs = ylbs = array("",Ncomb)
	envcex = speccex = 1 * sqrt(totDim-1); smpcex=1.5 * sqrt(totDim-1);
	
	PubQual = opt$report$PublicationQuality
	labelcex=1
	
	if (PubQual){
		mainPX = 800; 
		bty = "l"
	} else {
		mainPX = 800;
		bty = "o"
	}
	#stop("group in plot")
	CCAb=FALSE
	if (inherits(ef,"anova.cca") || inherits(mds,"cca")){
		CCAb=TRUE
		pv = ef[,4]; names(pv) = dimnames(ef)[[1]];
		smplDes = dimnames(mds$CA$u)[[1]]
		featDes = dimnames(mds$CA$v.eig)[[1]]
		sOb = summary(mds)
	} else {
		pv=c(ef$vectors$pvals,ef$factors$pvals);  
		names(pv)=c(names(ef$vectors$r),names(ef$factors$r));
		smplDes = dimnames(mds$points)[[1]]
		featDes = dimnames(mds$species)[[1]];
	}
	take = names(pv)[which(pv<=0.05)];
	hir=MetaInfo$hierarchy
	
	#Smp Descriptors
	
	if (is.null(colObj)){
		mSNC = getSampleColorName(smplDes,MetaInfo)
	} else {
		mSNC = colObj		
	}
	smpCol = mSNC$SampleColor;
	smpLab = mSNC$SampleNames;
	smpPch = mSNC$SamplePCH;
	colLeg = mSNC$ColLegend
	colGrad = mSNC$ColGradient
	#display only subset of samples?
	smplSel=array(TRUE,length(smplDes))
	if (!is.null(displaySubs)){
		smplSel=displaySubs[smplDes]
	}
	
	#label priority
	featPrio = array(1,length(featDes))
	if (!is.null(MetaInfo$FeatureImportance)){
		mm=match(featDes,names(MetaInfo$FeatureImportance))
		featPrio = order(MetaInfo$FeatureImportance[mm],decreasing =F)
	}
	
	#-------------------------------------------------
	#--- Start of big loop around all dimensions   ---
	if (totDim>2){labelcex=1.5}
	addiT = ""	
	if (CCAb){
		#variance in constrained axis
		eigs = c(sOb$concont[[1]][1,],sOb$cont[[1]][1,])
		#variance in constrained + unconstrained axis
		#eigs_all = sOb$cont[[1]][1,]
		#sOb$pcSDV = eigs
		addiT = c(array("constr.",length(sOb$concont[[1]][1,])),
				array(" unconstr.",length(sOb$cont[[1]][1,])))
		siteTag = "wa"
	} else {
		eigs=mds$eig; addiT = array("",length(eigs));
		siteTag = "sites"
	}
	for (dims in 1:Ncomb){
		sel1=  combis[dims,1]; sel2 = combis[dims,2]
		xlbs[dims] = paste(ordiMethod," D",sel1," (",formatC(eigs[sel1]/sum(eigs)*100),"%",addiT[sel1],")",sep="");
		ylbs[dims] = paste(ordiMethod," D",sel2," (",formatC(eigs[sel2]/sum(eigs)*100),"%",addiT[sel2],")",sep="");
	}
	#choices =combis[dims,] 
	
	if (IndFit){
		#-------------  Feature Plot  --------------------
		
		#colored hierarchical phylogeny plot
		if(!is.null(hir)){
			#browser()
			mm = match(dimnames(mds$species)[[1]],hir$names);
			for (j in 1:2){
				#hir[,1]=="Bacteria"
				#for (k in 1:j){
				curN = dimnames(hir$subcol)[[2]][j]
				writeInfo(myPaste("Feature Distribution - ",curN),Opt,"picf")
				plot(mds, dis = "sp", type = "n",xlab=xlbs,ylab=ylbs,main=curN)
				abline(h=0,v=0,lty=2,col="gray")
				spcol = hir$subcol[mm,j]
				plot.ef(ef, p.max = 0.8, col="black",xlab=xlbs,ylab=ylbs,noArrow=T,
						cex=envcex);
				points(mds,dis="sp",col=spcol,pch="+",cex=1)
				writeInfo("",Opt,"dev.off");
				writeInfo(myPaste("Legend - ",curN),Opt,"picf")
				plot(1)
				if(j>hir$subPos && length(hir$s.cols[[j]]) > 1){
					lcol = c(hir$s.cols[[j]][[1]],hir$s.cols[[j]][[2]]);
					llab = c(hir$s.lvls[[j]][[1]],hir$s.lvls[[j]][[2]]);
				} else {
					lcol = c(hir$s.cols[[j]][[1]]);
					llab = c(hir$s.lvls[[j]][[1]]);
					
				}
				legend("topright",llab,col=lcol,pch="+",ncol=ceiling(length(llab)/20))
				dev.off()
				#}
			}
		} 
		
		if (FeatureDistrPlot){
			#-----------------    Feature Distribution plot   --------------
			writeInfo("Feature Distribution",Opt,"picf",second=c(mainPX,Ncomb*mainPX))
			par(mfrow=c(1,Ncomb))
			for (dims in 1:Ncomb){
				plot(mds, dis = "sp", type = "n", bty=bty,
						xlab=xlbs[dims],ylab=ylbs[dims],cex=smpcex,choices =combis[dims,])
				abline(h=0,v=0,lty=2,col="gray")
				#browser()     orditorp.2
				sel <- orditorp(mds, dis = "sp", lab = FeatNam,col="black", pcol = "black", pch = "+",cex=speccex,priority=featPrio,choices =combis[dims,])
				#plot.ef(ef, p.max = 0.8, col="royalblue",xlab=xlbs,ylab=ylbs,noArrow=T);
				if (insertSitToSpecP){
					xtbl=scores(mds,display=siteTag,choices =combis[dims,],directed=FALSE)
					if(usePoints== -1){
						#text(mds, dis = siteTag, lab = smpLab,col=smpCol, cex=smpcex,choices =combis[dims,],xpd=NA)
						text(xtbl[smplSel,1],xtbl[smplSel,2], lab = smpLab[smplSel],col=smpCol[smplSel], cex=smpcex,xpd=NA)
					} else if(usePoints==0){
						sel <- orditorp(mds, labels=smpLab, dis = siteTag,  pcol = smpCol,col=smpCol, pch = "+",cex=smpcex,choices =combis[dims,],xpd=NA)
					} else { #publication quality: use symbols
						points(xtbl[smplSel,1],xtbl[smplSel,2],pch=smpPch[smplSel],col=smpCol[smplSel],bg=smpCol[smplSel],cex=1.5)
						if (is.null(colGrad)){
							legend("topleft",legend=names(colLeg),pch=attr(colLeg,"pch"),col=colLeg,pt.bg=colLeg)
						} else {
							legendGrad(colGrad)
						}
					}
					if (!is.null(groups)){
						addOrdiSampleGroups(xtbl[smplSel,],groups[smplSel],smpCol[smplSel])					
					}
				}
			}
			writeInfo("",Opt,"dev.off");
		}
		
		
		if (!is.null(ef) && (!CCAb || !is.na(scores(mds, c(1,2), "cn")	))){
			#----------  MetaData plot  ----------------
			writeInfo("metaData Correlations",Opt,"picf",second=c(mainPX,Ncomb*mainPX))
			par(mfrow=c(1,Ncomb))
			for (dims in 1:Ncomb){
				#store later in taken, who was selected, so the legend doesn't show empty categories
				sigs=c(0.01,0.05,0.1,0.5,1); taken = array(TRUE,length(sigs));
				colors2 = (rainbow(length(sigs),start=0,v=0.7)); colors = rainbow(length(sigs),alpha=0.25,v=0.7);
				main=paste("Metadata ",totDim,"D ",ordiMethod," fit",whichDim(combis,dims),sep="")
				if (CCAb){
					plot(mds,dis="cn",xlab=xlbs[dims], main=main,
							ylab=ylbs[dims])
					#pts = text(mds,dis="cn",select=c("GenoT_B6"))
				} else {
					#cex higher in this plot
					plot.ef(ef, noArrow=F,add=F,xlab=xlbs[dims], main=main,
							ylab=ylbs[dims],txtcol=colors2[1], arrow.mul=1,
							makeEmpty=T,cex=smpcex,choices =combis[dims,],
					);					
					
					#points(mds, dis = "sites", col="black",pch="x")
					for (i in 1:length(sigs)){
						if (i == 1) {ssig = -1} else {ssig = sigs[i-1]}
						sel  = names(pv)[pv<sigs[i] & pv > ssig]
						#print(sel)
						if (all(FALSE==sel)){
							taken[i] = FALSE; next;	}
						plot.ef(ef,takeI=sel, col=colors[i],noArrow=F,add=T,
								txtcol=colors2[i], arrow.mul=1,
								choices =combis[dims,],cex=envcex);
					}
					legend("topright", paste(formatC(paste("<",sigs[taken]))), cex=1, col=colors2[taken],pch=1)
					#points(mds)
					abline(h=0,v=0,lty=2,col="gray")
				}
			}
			writeInfo("",Opt,"dev.off");
			
		}
	}

	
	#temporarily deactivate surface fit:
	take = c()
	
	if (length(take)>3){take=take[  order(pv[pv<0.05])[1:min(maxTake,length(take))]];}
	
	writeInfo(paste(ordiMethod,"main Triplot - scaling 2"),Opt,"pic",second=c(mainPX,Ncomb*mainPX))
	par(mfrow=c(1,Ncomb))
	colors = array(1,length(class[[1]])); colors2 = colors-.5;
	
	mainDis = c("sp","si","cn") 
	if (justSites){mainDis="si"}
	for (dims in 1:Ncomb){
		wDim = whichDim(combis,dims);
		mainT = paste("Triplot - ",totDim,"D ",ordiMethod," -",wDim,sep="");
		p=plot(mds, dis=mainDis,type="n",main=mainT,xlab=xlbs[dims],bty=bty,
				ylab=ylbs[dims],cex=smpcex,choices =combis[dims,],
				cex.main=1.7,cex.axis=labelcex,cex.lab=labelcex)	
		xtbl=scores(mds,display=siteTag,choices =combis[dims,],directed=FALSE)
	
		abline(h=0,v=0,lty=2,col="gray")
		os = ordispace(externalFacHull,TimeSeries,hullType,mds,combis,dims,Opt,MetaInfo,
				colObj=mSNC,ordiSpCol=ordiSpCol,usePoints=usePoints,xtbl=xtbl,
				dispSubset= smplSel,showLbl=doXtraLegend)
		usePoints = os$usePoints
		spiderLgnd = NULL #only needs to be done once
		
			
		if (length(take) > 1 && IndFit){
			numClass = lapply(class,is.numeric);
			colors = rainbow(length(take));
			for (i in 1:length(take)){
				if (numClass[[take[i]]]){
					tmp <- ordisurf(mds, class[,take[i]], add = TRUE,col=colors[i],choices =combis[dims,]);
				}
			}
			legend("topleft", paste(take,formatC(pv[take])), cex=1, col=colors,pch=1)
		} else {
			if (insertMetaToSit){
				if (ordiMethod %in%c("dbRDA","RDA")){
					if (!is.null(mds$CCA$centroids)){
						text(mds,cex=envcex,choices =combis[dims,],dis="cn",xpd=NA)
					}
				} else {#NMDS etc
					plot.ef(ef, p.max = 0.8, col="darkgray",noArrow=T,add=T,
							xlab=xlbs[dims],ylab=ylbs[dims],choices =combis[dims,],
							cex=envcex);					
				}
			}
		}
		
		if(usePoints== -1){
			#text(mds, dis = siteTag, lab = smpLab,col=smpCol, cex=smpcex,choices =combis[dims,],xpd=NA)
			text(xtbl[smplSel,1],xtbl[smplSel,2], lab = smpLab[smplSel],col=smpCol[smplSel], cex=smpcex,xpd=NA)
		} else if(usePoints==0){
			sel <- orditorp(mds, labels=smpLab, dis = siteTag,  pcol = smpCol,col=smpCol, pch = "+",cex=smpcex,choices =combis[dims,],xpd=NA)
		} else if (usePoints==-2){
		} else { #publication quality: use symbols
			points(xtbl[smplSel,1],xtbl[smplSel,2],pch=smpPch[smplSel],col=smpCol[smplSel],cex=1.5)#bg=smpCol[smplSel],
			if (is.null(colGrad)){
				legend("bottomright",legend=names(colLeg),pch=attr(colLeg,"pch"),col=colLeg,pt.bg=colLeg)
			} else {
				legendGrad(colGrad)
			}
		}
		
		if (!is.null(groups)){
			addOrdiSampleGroups(xtbl[smplSel,],groups[smplSel],smpCol[smplSel])					
		}
		
		blowup = TRUE
		if (!justSites && (inherits(mds,"rda") || !is.na(mds$species))){
			xtbl2=scores(mds,display="sp",choices =combis[dims,],directed=FALSE)
			if (blowup){
				xtbl2[,1] = xtbl2[,1] *.5* diff(range(xtbl[,1])) / diff(range(xtbl2[,1])) 
				xtbl2[,2] = xtbl2[,2] *.5* diff(range(xtbl[,2])) / diff(range(xtbl2[,2])) 
			}
			text(xtbl2, labels=FeatNam,  col="darkred", cex=1)
		}
	}
	writeInfo("",Opt,"dev.off");
	
	return(invisible(NULL))
}
#'short funciton to create explanatory string for plot titles that are multidimnesional
whichDim = function(combis,dims){
	paste(" Dim ",combis[dims,1]," & ",combis[dims,2],sep=""); 	
}
#' NMDS helper: plots environmental fit object from vegan
#' @param x = ef object
plot.ef = function (x, choices = c(1, 2), arrow.mul, at = c(0, 0), axis = FALSE, 
		p.max = 1, col = "blue", add = TRUE, noArrow = FALSE, txtcol=NULL, 
		takeI = NULL ,makeEmpty = FALSE,cex=1, ...) 
{
	#browser()
	formals(arrows) <- c(formals(arrows), alist(... = ))
	if (is.null(txtcol)){txtcol = col;}
	vect <- NULL
	if (!is.null(p.max) || !is.null(takeI)) {
		if (!is.null(x$vectors)) {
			#browser()
			if (is.null(takeI)){	take <- x$vectors$pvals <= p.max
			} else { 
				subt=takeI %in%dimnames(x$vectors$arrows)[[1]];
				take = takeI[subt]; }
			
			x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
			x$vectors$r <- x$vectors$r[take]
			if (nrow(x$vectors$arrows) == 0) 
				x$vectors <- vect <- NULL
		}
		if (!is.null(x$factors)) {
			if (is.null(takeI)){	take <- x$factor$pvals <= p.max
			} else {
				matchTo = dimnames(x$factors$centroids)[[1]];
				take=grepl(paste("^",takeI[1],sep=""),matchTo)
				for (tmp in 2:length(takeI)){take = take | grepl(paste("^",takeI[tmp],sep=""),matchTo);}
				if (length(take) == 0){take=array(FALSE,length(x$factors$pvals))} 
			}
			x$factors$centroids <- x$factors$centroids[take, 
					, drop = FALSE]
			if (nrow(x$factors$centroids) == 0) 
				x$factors <- NULL
		}
	}
	if (!is.null(x$vectors)) {
		vect <- sqrt(x$vectors$r) * x$vectors$arrows[, choices, 
				drop = FALSE]
		if (missing(arrow.mul)) {
			if (!add) 
				arrow.mul <- 1
			else arrow.mul <- vegan:::ordiArrowMul(vect, at = at)
		}
		if (axis) {
			maxarr <- round(sqrt(max(x$vectors$r)), 1)
			ax <- -c(-1, 0, 1) * arrow.mul * maxarr
		}
		vect <- arrow.mul * vect
		vtext <- sweep(1.1 * vect, 2, at, "+")
		vect <- sweep(vect, 2, at, "+")
	}
	if (!add) {
		xlim <- range(at[1], vect[, 1], x$factors$centroids[, 
						1]) * 1.15
		ylim <- range(at[2], vect[, 2], x$factors$centroids[, 
						2]) * 1.15
		if (!is.null(vect)) 
			plot(vect, xlim = xlim, ylim = ylim, asp = 1, type = "n", 
					...)
		else if (!is.null(x$factors)) 
			plot(x$factors$centroids[, choices, drop = FALSE], 
					asp = 1, xlim = xlim, ylim = ylim, type = "n", 
					...)
		else stop("Nothing to plot")
	}
	if (makeEmpty){return(NULL)}

	if (!is.null(vect)) {
		if (!noArrow){
			arrows(at[1], at[2], vect[, 1], vect[, 2], len = 0.05, 
					col = col)
		}
		text(vtext, rownames(x$vectors$arrows), col = txtcol, cex=cex,...)
	}
	if (!is.null(x$factors)) {
		text(x$factors$centroids[, choices, drop = FALSE], rownames(x$factors$centroids), 
				col = txtcol, cex=cex, ...)
	}
	if (axis && !is.null(vect)) {
		axis(3, at = ax + at[1], labels = c(maxarr, 0, maxarr), 
				col = col)
		axis(4, at = ax + at[2], labels = c(maxarr, 0, maxarr), 
				col = col)
	}
	invisible()
}

addColMat.numeric = function(mat,lab,HQR,name,LOG=FALSE){
	#browser();
	ret = list();
	if (any( dimnames(mat)[[2]] == name) ){#already added some way or another
		ret$mat = mat; ret$lab = lab;
		return(ret);
	}
	matchTo = dimnames(mat)[[1]];
	nn = names(HQR);
	mm = match(matchTo,nn);
	HQR = HQR[mm];
	
	if(!LOG){
		HQRadd = as.numeric(HQR)
	} else {
		HQRadd = log10(as.numeric(HQR))
	}
	mat = cbind(mat,HQRadd);
	d2=dim(mat)[2];
	dimnames(mat)[[2]][ d2 ] = name
	HQR=getQuantsStr(HQR,c(.025,.5,.975))
	lab = cbind(lab,HQR);
	dimnames(lab)[[2]][ d2 ] = name
	
	ret$mat = mat; ret$lab = lab;
	return(ret);
}

auto.label.colRowCoLab = function(mat,lab)
{
	for (i in 1:dim(mat)[2]){
		if (all(lab[,i] == "") && all(!grepl("#[A-F0-9]{6}",mat[,i])) ){ #ok, time to find useful label
			uniqueLab = unique(mat[,i]);
			if (length(uniqueLab) < (dim(mat)[1]/3*2)) { #don't want complete cluter
				lab[,i] = mat[,i]
			}
		}
	}
	return(lab)
}

