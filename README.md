# Three ways to perform EOF analysis in MATLAB

Here we show three ways to perform EOF analysis using MATLAB codes.

## 0. Data
In this example, we use monthly 3D sst anomaly data during 1979-2020 to perform the EOF analysis. Data are stored as the `ssta` mat file, containing lon, lat, and 3D [91,31,504] sst anomaly data.

## 1. Direct calculation based on EOF essence
The essence of EOF analysis can be covered using following steps:
Consider a data `F`, where each row corresponds to a time point and each column corresponds to a spatial grid.

Step 1. Calculating covariance matrix

$C = cov(F)$

Step 2. Calculating eigenvalues and eigenvectors

$EOF,L=eig(C)$

Step 3. Calculating principal component time series

$PC=EOF^T \times (F-Mean)$

We can directly perform EOF analysis in MATLAB based on the EOF essence.

Before the application of EOF, the input SST data should be scaled by its latitude. 
```
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);
```
There, we need to reshape the 3D data (lon-lat-time) into 2D data (time-grid).
```
ssta=(reshape(ssta,91*31,504))';
```
The SST data contain `NaN`, corresponding to lands. We need to identify lands within grids and store them for further application.
```
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);
```
Now the SST data contain only SST data over oceans, without `NaN`.

Now we remove the mean in each time point of the SST data and calculate the covariance matrix.
```
F=detrend(ssta,0);
C=F'*F;
```
Then we calculate the eigenvectors and eigenvalues, which are EOFs and explained variance, respectively. The PCs are calculated based on matrix multiplication between EOFs and F.
```
[EOFs,D]=eig(C);
PCs=EOFs'*F';
```
After normalization of the eigenvalues D, they are confined to the range [0,1], which provides the proportion of explained variance in a percentage format.
```
D=diag(D);
D=D./nansum(D);
```
Then we identify the first EOF mode and PC time series, demonstrating the largest explained variance.
```
EOF1=EOFs(:,D==nanmax(D));
PC1=PCs(D==nanmax(D),:);
```
Then we reshape the first EOF mode back to the 2D format, using the land index `nanindex`. 
```
sEOF1=NaN(91*31,1);
sEOF1(~nanindex)=EOF1;
sEOF1=reshape(sEOF1,91,31);
```
After that, the EOF and PC should be scaled to standard units. 
```
sEOF1=sEOF1.*nanstd(PC1);
PC1=PC1./nanstd(PC1);
```
Then we write some codes to visualize the results.
```
figure
subplot(2,1,1);
m_proj('miller','lon',[nanmin(lona) nanmax(lona)],'lat',[nanmin(lata) nanmax(lata)]);
m_contourf(lona,lata,sEOF1',linspace(-1.4,1.4,200),'linestyle','none');
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('linewidth',2,'fontname','consolas');
colormap(m_colmap('diverging'));
caxis([-1.4 1.4]);
s=colorbar('fontname','consolas','fontsize',12);
title(s,'^{o}C','fontname','consolas');
set(gca,'fontsize',12)
title('EOF1: 34.92%','fontsize',16,'fontname','consolas');

subplot(2,1,2);
plot(1:504,PC1,'r','linewidth',2);
set(gca,'xtick',[6:60:504],'xticklabels',1980:5:2021,'fontname','consolas','fontsize',12);
xlabel('Year','fontname','consolas');
ylabel('PC1','fontname','consolas');
xlim([1 504]);
set(gca,'fontsize',12,'linewidth',2)
title('PC1: 34.92%','fontsize',16,'fontname','consolas');
```
![Image text](https://github.com/ZijieZhaoMMHW/Three_EOF/blob/main/EOFessence.png)
The first EOF mode shows an evidently positive ENSO phases, with a PC time series demonstrating significant correlation with the traditionally defined Nino34 index.

## 2. Calculation based on Principal Component Analysis (PCA)
The EOF is technically the same as the PCA. In fact, The term "spatial PCA" was initially used to refer to the EOF. Therefore, the EOF analysis can be achieved by using the MATLAB built-in function `pca` for principal component analysis.

Similar to the first step, we need to reconstruct the 3D data into 2D format. 
```
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);

ssta=(reshape(ssta,91*31,504))';
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);
```
Then we apply the `pca` function to the 2D SST data.
```
F=detrend(ssta,0);
[coef,score,latent]=pca(F);
latent=latent./nansum(latent);
```
The main outputs of `pca` include `coef`,`score`,and  `latent`. `coef` is the loadings from the PCA, which is just the EOF spatial pattern. `score` is the princpial component series projected to the `coef`. `latent` is principal component variance, which can be used as explained variance for each mode after normalization. 

Then we reshape the first EOF mode back to 2D format and visualize them after standardization.
```
scoef1=NaN(91*31,1);
scoef1(~nanindex)=coef(:,1);
scoef1=reshape(scoef1,91,31);

score1=score(:,1);
scoef1=scoef1.*nanstd(score1);
score1=score1./nanstd(score1);

figure
subplot(2,1,1);
m_proj('miller','lon',[nanmin(lona) nanmax(lona)],'lat',[nanmin(lata) nanmax(lata)]);
m_contourf(lona,lata,scoef1',linspace(-1.4,1.4,200),'linestyle','none');
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('linewidth',2,'fontname','consolas');
colormap(m_colmap('diverging'));
caxis([-1.4 1.4]);
s=colorbar('fontname','consolas','fontsize',12);
title(s,'^{o}C','fontname','consolas');
set(gca,'fontsize',12)
title('EOF1: 34.92%','fontsize',16,'fontname','consolas');

subplot(2,1,2);
plot(1:504,score1,'r','linewidth',2);
set(gca,'xtick',[6:60:504],'xticklabels',1980:5:2021,'fontname','consolas','fontsize',12);
xlabel('Year','fontname','consolas');
ylabel('PC1','fontname','consolas');
xlim([1 504]);
set(gca,'fontsize',12,'linewidth',2)
title('PC1: 34.92%','fontsize',16,'fontname','consolas');
```
![Image text](https://github.com/ZijieZhaoMMHW/Three_EOF/blob/main/EOFpca.png)

## 3. Calculation using the Climate Data Toolbox
The application of the EOF can also be achieved by using some peripherals, such as the [Climate Data Toolbox](https://chadagreene.com/CDT/CDT_Getting_Started.html) developed by [Chad Greene](https://github.com/chadagreene). In this toolbox, there is a function `eof` to directly performe EOF analysis on 3D spatiotemporal dataset. For this example, it can be achieved by following codes:
```
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);
[eof_maps,pc,expvar]=eof(ssta);
eof1=eof_maps(:,:,1);
pc1=(pc(1,:))';
eof1=eof1.*nanstd(pc1);
pc1=pc1./nanstd(pc1);

figure
subplot(2,1,1);
m_proj('miller','lon',[nanmin(lona) nanmax(lona)],'lat',[nanmin(lata) nanmax(lata)]);
m_contourf(lona,lata,eof1',linspace(-1.4,1.4,200),'linestyle','none');
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('linewidth',2,'fontname','consolas');
colormap(m_colmap('diverging'));
caxis([-1.4 1.4]);
s=colorbar('fontname','consolas','fontsize',12);
title(s,'^{o}C','fontname','consolas');
set(gca,'fontsize',12)
title('EOF1: 34.92%','fontsize',16,'fontname','consolas');

subplot(2,1,2);
plot(1:504,pc1,'r','linewidth',2);
set(gca,'xtick',[6:60:504],'xticklabels',1980:5:2021,'fontname','consolas','fontsize',12);
xlabel('Year','fontname','consolas');
ylabel('PC1','fontname','consolas');
xlim([1 504]);
set(gca,'fontsize',12,'linewidth',2)
title('PC1: 34.92%','fontsize',16,'fontname','consolas');
```
![Image text](https://github.com/ZijieZhaoMMHW/Three_EOF/blob/main/EOFcdt.png)
