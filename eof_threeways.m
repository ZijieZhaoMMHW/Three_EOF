%% Reformat Data
lon(lon<0)=360+(lon(lon<0));
[lon_re,i]=sort(lon);

sst_re_m=sst_anom_m(i,:,:);

ssta=sst_anom_m(lon_re<=300 & lon_re>=120,lat<=30 & lat>=-30,:);
lona=lon_re(lon_re<=300 & lon_re>=120);
lata=lat(lat<=30 & lat>=-30);

%% EOF - Direct Calculation
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);

ssta=(reshape(ssta,91*31,504))';
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);

F=detrend(ssta,0);
C=F'*F;

[EOFs,D]=eig(C);
PCs=EOFs'*F';

EOF1=EOFs(:,end);
PC1=PCs(end,:);

sEOF1=NaN(91*31,1);
sEOF1(~nanindex)=EOF1;
sEOF1=reshape(sEOF1,91,31);

sEOF1=sEOF1.*nanstd(PC1);
PC1=PC1./nanstd(PC1);

%% EOF - CDT
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);
[eof_maps,pc,expvar]=eof(ssta);
eof1=eof_maps(:,:,1);
pc1=(pc(1,:))';
eof1=eof1.*nanstd(pc1);
pc1=pc1./nanstd(pc1);

%% EOF - PCA
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);

ssta=(reshape(ssta,91*31,504))';
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);

F=detrend(ssta,0);
[coef,score,latent]=pca(F);

scoef1=NaN(91*31,1);
scoef1(~nanindex)=coef(:,1);
scoef1=reshape(scoef1,91,31);

score1=score(:,1);
scoef1=scoef1.*nanstd(score1);
score1=score1./nanstd(score1);
