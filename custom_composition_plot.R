#' species/Cog/Feature composition graph
#' @param labeling: which kind of labeling to make underneath plot c(colorBar,Legend)
#' @param order.by: either a vector with specific orderings, a cluster object (to project a cluster next to composition) or a keyword ("SampleCol",NULL)
#' @param as.percentage: display read counts or % of total composition on x-axis?
#' @param disp: plot type (polygon,bar)
#' @param SampleGroupLink: links two groups that should be displayed together
#' @param grPix set pixels of graph
Interface.SpeciesComposition = function(Mat,Opt=extractOptions(Mat),min.include=10,cexLabRow=1,
		title="Composition of Samples", labeling = "Legend",matTy="Normed",
		CatName="",totcex=min(cexLabRow/1.2,1),grPix=NULL ,
		MetaInfo = extractMetaInfo(Mat), overlayMat= NULL, 
		plotImportance=FALSE,labelSpace=0.4,
		FeatureSubset=NULL, poolFeatures=NULL, subset=NULL,as.percentage=FALSE,
		order.by=NULL,labCol=NULL,higherLvl=NULL,disp=c("bar","polygon"),
		TimeSeries=NULL,SampleGroups=NULL,SampleGroupLink=NULL){
	#browser()
	disp = match.arg(disp)
	M_ = getMatrix(Mat,type="Normed",row.subset=FeatureSubset,col.subset=subset,poolFeatures=poolFeatures,getRest=TRUE,
			filtered=F,higherLvl=higherLvl)
	matWarn(M_,matIs="Trans");
	
	unit = ""
	if(as.percentage){
		M_=t(t(M_)/colSums(M_)) * 100
		unit = "%"
	}
	if (!is.null(TimeSeries)){
		disp="mix"
		if (!is.null(order.by)){
			writeInfo("Sample Order ignored as TimeSeries argument is provided",Opt,"warning")
		}
		if (is.null(SampleGroups)){
			order.by=order(TimeSeries)		
		} else {#include sample groups info in ordering
			#match(names(TimeSeries),names(SampleGroups))
			order.by=groupWiseTimeOrder(TimeSeries,SampleGroups,TRUE,SampleGroupLink,TRUE)
			groupOrdering = attr(order.by,"groupOrder")
		}
	}
	
	# sample ordering: cluster preparation
	SampleClus = NULL;
	if (!is.null(order.by) && inherits(order.by,"DDT.unsuper.cluster")){
		#myLibs("gclus")
		#otree = gclus::reorder.hclust(otree,dd)
		if (inherits(order.by$clusterObj,"hclust")){
			ord = order.by$clusterObj$order
			names(ord) = order.by$clusterObj$labels
			SampleClus = order.by$clusterOb
			order.by = ord
		} else {
			order.by = names(order.by$cluster)[order(order.by$cluster)]
		}
	}
	
	#fine tuning of graph parameters
	#if (!is.null(TimeSeries) && length(TimeSeries)>50){
		#labelSpace = .2; #.2
	#} else 	if (Opt$report$PublicationQuality){
		#labelSpace = .4; #.2
		#pixels = c(700,600) #700,600 c(opt$report$globalRes,opt$report$globalRes)
	#}		#browser()
	
	pixels = c((dim(M_)[2]*25+40)*(1+labelSpace),600)
	if (!is.null(grPix)){
		pixels=grPix
	}
	
	#first cut matrix to X most abundant species
	min.include = min(dim(M_)[1], min.include)
	meanSp = rowMeans(M_);
	ordMS = order(meanSp,decreasing=TRUE);
	if (!is.null(overlayMat)){ ovM = overlayMat[ordMS[1:min.include],]; }
	M = M_[ordMS[1:min.include],];
	if (dim(M)[1]!=dim(M_)[1]){
		M=rbind(M,apply(M_[ordMS[-(1:min.include)],,drop=F],2,sum)   )
		dimnames(M)[[1]][min.include+1] = "other"
		if (!is.null(overlayMat)){
			ovM=rbind(ovM,apply(overlayMat[ordMS[-(1:min.include)],],2,sum)   )
			dimnames(ovM)[[1]][min.include+1] = "other"
		}
	}
	othmatch = which("other"==dimnames(M)[[1]])
	if (length(othmatch) >= 1){
		othAbund = colSums(M[othmatch,,drop=FALSE])
		M = M[-othmatch,,drop=FALSE]
		M = rbind(M,other=othAbund)
	}
	
	# set colors
	nc = ncol(M); nr = nrow(M);
	
	mSNC = getSampleColorName(M,MetaInfo,givenCol=labCol);
	labRowCol = mSNC$SampleColor;
	labRow = mSNC$SampleNames;

	#sort abundance data
	if (!is.null(order.by)){
		if (is.numeric(order.by) && length(order.by) == dim(M)[2]){
			if (!is.null(names(order.by))){
				mm = match (dimnames(M)[[2]],names(order.by))
			} else {
				mm = 1:length(order.by)
			}
			ordM = order.by[mm];
		} else if (length(order.by)==1 && order.by=="SampleCol"){
			lev = mSNC$ColLegend
			comM = matrix(0,dim(M)[1],length(lev))	
			for (i in 1:length(lev)){
				comM[,i] = rowMeans(M[,labRowCol==lev[i],drop=F])
			}
			#dimnames(comM)[[2]] = lev
			ord1 = order(comM[1,])
			lev1 = lev[ord1]
			ordh1 = array(NA,length(labRow))
			for (i in 1:length(lev1)){
				ordh1[ which(labRowCol==lev1[i]) ] = i
			}
			ordh2 = order(M[1,])
			ordM = order(ordh1,M[1,])
			
		} else if (length(order.by)==1){#chracter: match this
			sel = .extractComplexFeature(order.by,dimnames(M)[[1]])[[1]];
			if (sum(sel)==0){
				vals = M[1,]
			}
			vals = apply(M[sel,,drop=F],2,sum)
			ordM = order(vals)
		} else if (is.character(order.by)){
			ordM=order.by
			ordM =ordM[ ordM %in% dimnames(M)[[2]] ]
		}
		
	} else {
		ordM = order(M[1,])	
	}
	
	if (is.character(order.by)){
		namesOrd = ordM
	} else {
		namesOrd = dimnames(M)[[2]][ordM]
	}
	breakAfter=array(FALSE,dim(M_)[2])
	if (!is.null(SampleGroups)){
		mm=match(dimnames(M_)[[2]],names(SampleGroups))
		SampleGroups=SampleGroups[mm]
		if (all(SampleGroups==SampleGroups[1])){
			SampleGroups = NULL
		} else {
			breakAfter=SampleGroups[ordM]!=SampleGroups[ordM][c(2:length(mm),1)]
			#breakAfter [which(SampleGroups[order.by]!=SampleGroups[order.by][c(2:length(mm),1)]) ] = TRUE
			names(breakAfter) = dimnames(M_)[[2]]
		}
	}
	
	
	
	
	
	#---------- PLOT begins here --------------

	writeInfo("Species Composition",Opt,"pic",second=pixels,specialTreatment=TRUE)	
	ColLabWidth =  1.05*quantile(strwidth(labRow, cex = cexLabRow,"inc"),0.95); #I need at least 0.95 to fit in
	marin = c(ColLabWidth + ifelse(!is.null(SampleClus),0,0.25), .1);
	sidelabel=NULL
	#check if to plot SampleImportance
	plotSide = FALSE
	if (labeling == "Legend"){labelSpace = labelSpace*2;}
	lhei = c(4,labelSpace); lwid=c(4); lmat = matrix(0,2,1); 
	lmat[1,1]=1; lmat[2,1]=2; 
	if (!is.null(MetaInfo$SampleImportance) && plotImportance){
		plotSide = TRUE
		lwid = c(lwid,.17); lmat = cbind(lmat+1,NA); dlm = dim(lmat)[2];
		lmat[,dlm] = lmat[,dlm-1]
		lmat[lhei==4,dim(lmat)[2]] = 1
		
#		if (!is.null(MetaInfo$SampleImportance)){
			ColLabWidth =  1.05*quantile(strwidth(MetaInfo$SampleImportance, cex = cexLabRow,"inc"),0.95); #I need at least 0.95 to fit in
			marin[2] = ColLabWidth + 0.1;
			HQR = MetaInfo$SampleImportance;
			mm = match(dimnames(M)[[2]],names(HQR))
			HQR = HQR[mm];
			sideCol = colVec(HQR);
			sidelabel = HQR;
			locsrt = 0;locadj=c(0.5,1)
			locsrt = 0;locadj=c(0,.5)
#		}
		
	}
	if(!is.null(title)){
		lhei = c(.2,lhei); lmat = rbind(1,lmat+1)
	}
	if (!is.null(SampleClus)){
		lwid = c(.3,lwid); lmat = cbind(max(lmat)+1,lmat);
		lmat[2,1]= lmat[2,2]
	}
	lmat[is.na(lmat)] <- 0
	nlevs = dim(M)[1]
	#getGeneralColors(nlevs,opt)
	usedCols <- getGeneralColors(nlevs,Opt)
	#usedCols[c(-1,-4)] = grayify(usedCols[c(-1,-4)],2)
	x=layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	
	#title
	if(!is.null(title)){
		par(mar=c(0,0,0,0));
		plot.new();
		#limit locex to min of 5, character height or char width (relative to available space)
		locex = min(5,  min(par("fin")[1]/ strwidth(title,units="inch") ,par("fin")[2]/ strheight(title,units="inch"))) 
		if (!is.null(Opt$what)){
			title = paste(title," (",Opt$what,")",sep="");
		}
		text(.5,.4,adj=c(0.5,0.5),lab=title,cex=locex*.8,xpd=NA)
	}
	
	M = M[,ordM]; if (!is.null(overlayMat)){ ovM = ovM[,ordM];}
	if (dim(M)[2]<20){
		space = 0.1
	} else {
		space = 0
	}
	#side color
	if (plotSide && plotImportance) {
		HQR = MetaInfo$SampleImportance;
		mm = match(dimnames(M)[[2]],names(HQR))
		HQR = HQR[mm];
		sideCol = colVec(HQR);
		sidelabel = array("",length(sideCol));
		sidelabel[which.min(HQR)] = min(HQR)
		sidelabel[which.max(HQR)] = max(HQR)
		mHQR = median(HQR);
		sidelabel[HQR==mHQR] = mHQR
		
		
		par(mai = c( 0, 0, 0, 0.1))
		plot.side.color(color = sideCol,label = sidelabel,horizontal=FALSE,
				SideLab = "Feature Sum")#, spaceBetween = space)
	}
	#real abundance plot

	par(mai=c(0,marin[1],0,marin[2]))
	#order for maximum abundance of first entry
	#ylims = c(0, ( (1+space)*dim(M)[1] -1) )
	
	
	if (disp%in%c("polygon","mix")){
		bp=polygonplot(M,col=usedCols,horiz=TRUE,breakA = breakAfter)
	} else {
		par(xaxs="r");	par(yaxs="r")		
		bp = barplot(M,col=usedCols,space=space,horiz=TRUE,border=NA,
			names.arg=array(NA,dim(M)[2]),axes = F,main=NA,sub=NA)#,ylim=c(1,11.1))
	}
	
	
	
	
	#get dimensions to plot tree later
	mainusr = par("usr")
	perc = grconvertY(min(bp),from="user",to="nfc")
	perc[2] = grconvertY(max(bp),from="user",to="nfc")
	incHei = par("fin")[2]
	
	
	##  --------------   insert overlay matrix
	if (!is.null(overlayMat)){
		xplo = yplo = matrix(0,prod(dim(M)),2)
		#colv = array("black",prod(dim(M)))
		dM1 = dim(M)[1];
		for (i in 1:dim(M)[2]){
			for (j in 1:dim(M)[1]){
				if (j==1){pxplo=0;} else {pxplo = sum(M[1:(j-1),i]);}
				idx = j+(i-1)*dM1
				xplo[idx,] = c(pxplo, pxplo + ovM[j,i])
				yplo[idx,] = c(bp[i],bp[i]);
				#colv[idx] = antiUsedCols[j]
			}
		}
		segments(x0=xplo[,1],y0=yplo[,1],x1=xplo[,2],y1=yplo[,2]
				,col="black",lty=1,lwd=1)
	}
	if (!is.null(sidelabel)){
		text(par("usr")[2],1:nc-(.5-space)+(1:nc*space), 
				srt=locsrt,labels = paste("",sidelabel[ordM]),
				col = sideCol[ordM],adj=locadj,xpd=TRUE,cex = cexLabRow)
	}
	
	#axis(2, 1:nc-.5, labels = labRow, las = 1, line = -0.7, tick = 0,cex.axis = cexLabRow,cex=totcex,col = labRowCol)
	#txtat = 1:nc-(.5-space)+(1:nc*space);
	txtat = bp
	subset = 1:(length(txtat))
	if (FALSE && !is.null(SampleGroups)){
		#attr(bp,"block_location")
		text(0,attr(bp,"block_location"),labels=paste(groupOrdering," "),adj=c(1,.5),xpd=TRUE,cex = cexLabRow)
	} else {
		text(0,txtat[subset], labels = paste(labRow[ordM][subset]," "),col = labRowCol[ordM][subset],adj=c(1,.5),xpd=TRUE,cex = cexLabRow)
	}
	tickAt = axTicks(1);
	shi = strheight("X","inch")
	lticks=length(tickAt)
	#browser()
	maxx =par("usr")[[2]]
	if (round(tickAt[lticks]) >= round(maxx-(maxx/25))){
		abline(v=tickAt[c(-1,-lticks)],lty=2,col="white",lwd=2)
	} else {	
		abline(v=tickAt[-1],lty=2,col="white",lwd=2)
	}
	subset = 1:(lticks)-1
	text(tickAt[subset],0,lab=paste(formatC(tickAt[subset]),unit,sep=""),xpd=NA,adj=c(0.5,1),cex=totcex);
	subset = (length(tickAt))
	text(tickAt[subset],0,lab=paste(formatC(tickAt[subset]),unit,sep=""),xpd=NA,adj=c(1,1),cex=totcex);
	
	#Color assignments
	mymai = c(0.1,.1,.02,.1)
	par(mai=mymai)
	#par(mar=c(1,2,2,1))
	#the "perfect" string assortment should at maximum take up 75% of horizontal space
#	lt = getFeatureNames(dimnames(M)[[1]],metaInfo)
	if (!is.null(poolFeatures)){
		lt = dimnames(M)[[1]]
	} else {
		lt = shortenHierNames(dimnames(M)[[1]],readable=Opt$report$PublicationQuality);	
	}
		
	
	
	# colored legend
	.legend.bar(lt,usedCols,CatName,antiColor=T,labeling=labeling,
			locex=1.2,strHeiInc = shi);
	
	#cluster on side?
	if (!is.null(SampleClus)){
		par(mai=c(0,0,0,0))
		loL = incHei * perc[1];
		upL = incHei * (1-perc[2]);
		
		par(mai=c(loL,.1,upL,0))
		plot(as.dendrogram(SampleClus), horiz = TRUE, 
				axes = FALSE, yaxs = "i", leaflab = "none"
				,ylim=c(1,length(SampleClus$labels)))
				#,ylim=(range(bp)+c(-.5,-.5)) )
	}
	
	
	
	writeInfo("",Opt,"dev.off");
	
	return(invisible(namesOrd))
	
}
