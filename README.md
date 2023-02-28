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

```


