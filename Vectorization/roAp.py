# Importing packages as needed
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astroquery.mast import Catalogs
from astroquery.utils.tap.core import Tap
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from sympy import *
from astropy.table import join, Table
import time
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import EllipseSelector, SpanSelector
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from sympy import *
import lightkurve as lk
import sys
import re
# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

class GetDataPlot :
    
    def __init__(self, objectName, objectRadius) :
        """
        Initializing the class with the object name and radius.
        """
        self.objectName = objectName
        self.radiusDeg = objectRadius

    def crossReferenceGaiaStarhorse(self, objectName, radius) :
        """
        Does a cone search query using Gaia and Starhorse,
        then cross reference the two and merge them. 
        """
        simbadResultTable = Simbad.query_object(objectName)
        
        if simbadResultTable is not None and 'ra' in simbadResultTable.colnames \
            and 'dec' in simbadResultTable.colnames:
            # Get the coordinates of the object
            simbadRA = simbadResultTable['ra'][0]
            simbadDEC = simbadResultTable['dec'][0]

            simbadCoord = coord.SkyCoord(simbadRA, simbadDEC, unit=(u.deg, u.deg), frame="icrs")

            # Define your query
            gaiaQuery = f"""
            SELECT 
                SOURCE_ID, ra, dec, parallax, pmra, pmdec,
                phot_g_mean_flux, phot_g_mean_flux_error, phot_g_mean_mag,
                phot_bp_mean_flux, phot_bp_mean_flux_error, phot_bp_mean_mag,
                phot_rp_mean_flux, phot_rp_mean_flux_error, phot_rp_mean_mag,
                bp_rp, bp_g, g_rp, l, b
            FROM
                gaiadr3.gaia_source
            WHERE
                CONTAINS(POINT('ICRS', gaiadr3.gaia_source.ra, gaiadr3.gaia_source.dec), 
                CIRCLE('ICRS', {simbadCoord.ra.deg}, {simbadCoord.dec.deg}, {radius})) = 1
            """
            gaiaJob = Gaia.launch_job_async(gaiaQuery)
            gaiaResult = gaiaJob.get_results()
            if len(gaiaResult) > 0:
                # Perform the crossmatch using TAP (Table Access Protocol)
                tapServiceURL = "http://tapvizier.u-strasbg.fr/TAPVizieR/tap/"
                tapService = Tap(url=tapServiceURL)
                # Selecting Starhorse catalog
                query = f"SELECT * FROM \"I/354/starhorse2021\" WHERE CONTAINS(POINT('ICRS', \"RA_ICRS\", \"DE_ICRS\"), \
                        CIRCLE('ICRS', {simbadCoord.ra.deg}, {simbadCoord.dec.deg}, {radius})) = 1"
                maxRetries = 3
                retryDelay = 10  # seconds
                for attempt in range(1, maxRetries + 1):
                    try:
                        starhorse_job = tapService.launch_job_async(query)
                        starhorseResult = starhorse_job.get_results()
                        break  # Break out of the loop if successful
                    except TimeoutError as e:
                        print(f"TimeoutError: {e}. Retrying after {retryDelay} seconds...")
                        time.sleep(retryDelay)
                else:
                    raise RuntimeError(f"Failed after {maxRetries} attempts. Aborting.")
            
                return gaiaResult, starhorseResult
            else:
                print("No Gaia sources found in the specified radius.")
                return None, None
        else:
            print(f"Simbad query for object {objectName} returned no results. Returning None.")
            return None, None

    def save_to_csv(self, results, filename):
        """
        Saves data to a csv file.
        """
        if results is not None:
            #df = results.to_pandas()
            results.to_csv(filename, index=False)
        else:
            print("Results are empty. Not saving to CSV.")
            
    def pmSelector(self, filename="Data/merged.csv", newpath='Data/pm-selected.csv'):
        """
        Graphs the data and lets you select a region
        to be filtered out, ridding of stars outside
        the region.
        """
        # Read in the merged data
        mergedData = pd.read_csv(filename)
        pmra = mergedData['pmra']
        pmdec = mergedData['pmdec']
        g = mergedData
        # Plot the proper motion diagram
        fig, ax = plt.subplots()
        ax.scatter(pmra, pmdec, s=1)
        plt.title("Proper Motion Density Plot")
        ax.set_xlabel("pmra")
        ax.set_ylabel("pmdec")
        
        def onselect(eclick, erelease):
            """
            Selects a region on the graph.
            """
            global xmin, xmax, ymin, ymax
            # Update the scatter plot with the selected region
            xmin = min(eclick.xdata, erelease.xdata)
            xmax = max(eclick.xdata, erelease.xdata)
            ymin = min(eclick.ydata, erelease.ydata)
            ymax = max(eclick.ydata, erelease.ydata)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            print(f'\nxmin:{xmin},xmax:{xmax},\nymin:{ymin},ymax:{ymax}')
            plt.draw()
                
        # Create an EllipseSelector
        selector = EllipseSelector(ax, onselect, interactive=True)
        
        plt.show()
        
        # Apply the condition to filter rows after the plot is closed
        if xmin is not None and xmax is not None and ymin is not None and ymax is not None:
            condition = ((pmdec > ymin) & (pmdec < ymax)) & ((pmra > xmin) & (pmra < xmax))
            df = g.loc[condition]

            # Save the DataFrame back to the CSV file
            df.to_csv(newpath, index=False)
            print(f'Selected data saved to {newpath} of length {len(df)} vs original length {len(g)}')
        else:
            print('No region selected.')
            
    def parallaxSelector(self, filename='Data/pm-selected.csv', newpath='Data/parallax-selected.csv'):
        """
        Graphs the data and lets you select a region
        to filter background stars out.
        """
        # Read in the merged data
        pmSelected = pd.read_csv(filename)
        p = pmSelected
        parallax = pmSelected['parallax']
        
        # Plot the proper motion diagram
        fig, ax = plt.subplots()
        ax.hist(parallax, bins='auto')
        plt.title("Parallax Histogram")
        ax.set_xlabel("parallax")
        ax.set_ylabel("count")
        
        def onselect(xmin, xmax):
            """
            Selects a region on the graph.
            """
            global pxmin, pxmax  # Declare the variables as global to modify them
            pxmin = xmin
            pxmax = xmax
            print(f'xmin: {pxmin}, xmax: {pxmax}')
            plt.draw()
                
        # Create an EllipseSelector
        selector = SpanSelector(ax, onselect, 'horizontal', useblit=True, interactive=True)
        
        plt.show()
        
        # Apply the condition to filter rows after the plot is closed
        if pxmin is not None and pxmax is not None:
            condition = (parallax > pxmin) & (parallax < pxmax)
            df = p.loc[condition]

            # Save the DataFrame back to the CSV file
            df.to_csv(newpath, index=False)
            print(f'Selected data saved to {newpath} of length {len(df)} vs original length {len(p)}')
        else:
            print('No region selected.')
        
        # Perform the parallax selection
    def importData(self):
        """
        Importing the known roAp csv file and
        the csv file from the 'GetData' class.
        """
        r = pd.read_csv('Data/field_roAp_gaiaNstarhorse.csv', index_col=False)
        c = pd.read_csv('Data/parallax-selected.csv', index_col=False)
        # Dataframe of each dataset
        self.roAp = pd.DataFrame(r)
        self.clust = pd.DataFrame(c)
        # Extracting features of all datasets
        self.yr = (self.roAp['GMAG0']).values
        self.yc = (self.clust['GMAG0']).values
        self.xr = (self.roAp['BP-RP0']).values
        self.xc = (self.clust['BP-RP0']).values
        
    def calculatingError(self):
        """
        Transposing the datasets and calculating the
        error ellipses using magnitude errors. Need
        to make this look nicer and break up into
        seperate functions!!
        """
        # For roAp dataset
        xr_array = np.array([self.xr])
        xrp_array = xr_array.transpose()
        yr_array = np.array([self.yr])
        yrp_array = yr_array.transpose()
        # For cluster dataset
        xc_array = np.array([self.xc])
        yc_array = np.array([self.yc])
        # Only transposing one so that we can evenly distribute the values
        xrc = (xrp_array-xc_array)
        yrc = (yrp_array-yc_array)
        # Finding the angle using arctan(opposite/adjacent)
        theta0 = np.arctan(yrc/xrc)
        theta0 = np.array(theta0)
        yd = yrc*(np.sin(theta0))
        xd = xrc*(np.cos(theta0))
        self.rd = (xd**2+yd**2)**(1/2) # Resultant vector for the distance between cluster and roAp stars
        
        # Error values for roAp and cluster
        gmre = np.array([(self.roAp['phot_g_mean_mag_error'])])
        bp_rp_re = np.array([((self.roAp['phot_bp_mean_mag_error']))-((self.roAp['phot_rp_mean_mag_error']))])
        self.gmce = np.array([(self.clust['phot_g_mean_mag_error'])])
        self.bp_rp_ce = np.array([((self.clust['phot_bp_mean_mag_error']))-((self.clust['phot_rp_mean_mag_error']))])
        # Squaring all terms
        sq_gmre = gmre**2
        sq_bp_rp_re = bp_rp_re**2
        sq_gmce = self.gmce**2
        sq_bp_rp_ce = self.bp_rp_ce**2
        # Transposing roAp array 
        gmre_p = (sq_gmre.transpose())
        bp_rp_rep = (sq_bp_rp_re.transpose())
        # Sigma squared terms
        gme = (gmre_p+sq_gmce)
        bp_rp_e = (bp_rp_rep+sq_bp_rp_ce)
        # Calculating sqrt(Sigma) from all errors
        sig_gme = (gme)**(1/2)
        sig_bp_rp_e = (bp_rp_e)**(1/2)
        # Making the Ellipse Equation
        a1 = sig_gme
        b1 = sig_bp_rp_e
        theta1 = np.arctan(a1/b1)
        theta1 = np.array(theta1)
        yD = a1*(np.sin(theta1))
        xD = b1*(np.cos(theta1))
        self.rD = (xD**2+yD**2)**(1/2) # Resultant vector for the edge of the ellipse
        
    def plotCMD(self):
        """
        Plotting the CMD with cluster
        and roAp stars.
        """
        fig, ax = plt.subplots(1,figsize=(6,6))
        # Plotting the cluster and roAp stars
        ax.scatter(self.xc,self.yc,c='m',label='Cluster',s=1, zorder=1)
        ax.scatter(self.xr,self.yr,c='c',label='roAp',s=2, zorder=2)
        
        plt.gca().invert_yaxis()

        plt.title(f'{self.objectName} and Known roAp Stars CMD')
        plt.ylabel('Absolute Mag [GMAG]')
        plt.xlabel('Color/Temp [Bp-Rp]')
        plt.grid(zorder=0)
        plt.legend(loc='best')
        plt.savefig(f'Data/CMD of {self.objectName}')
        plt.show()

    def getCandidates(self):
        """
        Using errors to find candidates in the
        cluster which was used for input, saving
        to a csv file.
        """
        self.plotCMD()
        
        gmcet = self.gmce.transpose()
        bp_rp_cet = self.bp_rp_ce.transpose()
        # R squared values for Cluster and roAp mag error
        self.cands = []
        addedDesignation = set()
        self.sigma = input('Amount of error to use : ') # User change the number of sigmas
        self.sigma = float(self.sigma)
        # Doing vectorization to classify candidates
        for q in range(len(self.xc)):
            for i in range(len(self.xr)):
                if (self.sigma*abs(self.rD[i][q])) >= (abs(self.rd[i][q])):
                    sourceID = int(self.clust.iloc[q]['SOURCE_ID'])
                    if sourceID not in addedDesignation: 
                        self.cands.append({'Designation':(self.clust['SOURCE_ID'])[q],'RA':(self.clust['ra'])[q],
                                        'DEC':(self.clust['dec'])[q],'GMAG':(self.clust['GMAG0'])[q],'GMAG_Error':(gmcet)[q],
                                        'GMAG_Sigma_Error':(self.sigma*gmcet)[q],'BP-RP':(self.clust['BP-RP0'])[q],
                                        'BP-RP_Error':(bp_rp_cet)[q],'BP-RP_Sigma_Error':(self.sigma*bp_rp_cet)[q],
                                        'Ratio-to-Sigma':((self.rd)[i][q])/((self.rD)[i][q])})
                        addedDesignation.add(sourceID)
                    else:
                        continue
        if self.cands == []:
            print('No candidates found.')
            sys.exit()
        else:
            print(f'Found {len(self.cands)} candidates.')
            # Saving the candidates to a csv file
            self.cands = pd.DataFrame(self.cands)
            self.cands.to_csv(f'Data/Candidates {self.objectName}.txt')
            self.cands.to_csv(f'Data/Candidates {self.objectName}.csv')
            self.cands.head(n=len(self.cands))
    
    def plottingCandidates(self):
        """
        Plotting the candidates with the
        roAp stars from the given cluster .
        """
        ub = self.cands['BP-RP']
        ua = self.cands['GMAG']
        xe = self.cands['BP-RP_Error']
        ye = self.cands['GMAG_Error']

        fig, ax = plt.subplots(1,figsize=(6,6))

        ellipse0 = [mpatches.Ellipse(xy=((ub)[p],(ua)[p]),width=(self.sigma*xe)[p],height=(self.sigma*ye)[p])
                for p in range(len(ua))]
        for e1 in ellipse0:
            ax.add_artist(e1)
            e1.set_facecolor(color='None')
            e1.set_edgecolor('#ff0000')
        e1.set_label(f'{self.sigma}$\sigma$ Cand Error')

        #plt.errorbar(ub,ua,yerr=c='red',label='Candidate',s=5,zorder=3)
        plt.scatter(x=ub,y=ua,c='m',label='candidates',s=3,zorder=3)
        plt.scatter(x=self.xr,y=self.yr,c='c',label='roAp',s=3,zorder=2)

        plt.gca().invert_yaxis()

        plt.title(f'Blue Stragler Candidates from {self.objectName}')
        plt.ylabel('Absolute Mag [GMAG]')
        plt.xlabel('Color/Temp [Bp-Rp]')
        plt.grid(zorder=1)
        plt.legend(loc='best')
        plt.savefig(f'Data/CMD with Cands of {self.objectName}')
        plt.show()
        
    def main(self):
        """
        Main function does the input for the query,
        cluster name and radius (in deg), and calculates the magnitude error
        used in calculations, then merges the errors to the cross-referenced data. 
        Then it calls the plotting functions to perform selection.
        Also, calls the function to find and plot the candidates.
        """
        self.objectName = self.objectName
        self.radiusDeg = self.radiusDeg
        # Perform cross-referencing with cone search
        gaiaResult, starhorseResult = self.crossReferenceGaiaStarhorse(self.objectName, self.radiusDeg)

        # Check if results are not empty
        if gaiaResult is not None and starhorseResult is not None:

            # Print column names for debugging
            print("Gaia Columns:", gaiaResult.colnames)
            print("StarHorse Columns:", starhorseResult.colnames)

            # Merge based on column names
            starhorseResult.rename_column('Source', 'SOURCE_ID')
            mergedResults = join(gaiaResult, starhorseResult, keys='SOURCE_ID', join_type='inner')
            
            # Add magnitude errors to Gaia data
            mergedResults['phot_g_mean_mag_error']  = (2.5 / np.log(10)) * (mergedResults['phot_g_mean_flux_error'] / mergedResults['phot_g_mean_flux'])
            mergedResults['phot_bp_mean_mag_error'] = (2.5 / np.log(10)) * (mergedResults['phot_bp_mean_flux_error'] / mergedResults['phot_bp_mean_flux'])
            mergedResults['phot_rp_mean_mag_error'] = (2.5 / np.log(10)) * (mergedResults['phot_rp_mean_flux_error'] / mergedResults['phot_rp_mean_flux'])

            # Save the merged results to a single CSV file
            mergedResults = mergedResults.to_pandas()
            #n_merged_results = merged_results.drop(['solution_id'], axis=1)
            self.save_to_csv(mergedResults, 'Data/merged.csv')

            # Display the merged results
            print("Merged Results:")
            print(mergedResults)

        else:
            sys.exit()
        
        # Perform the proper motion and parallax selection, finding candidates and plotting them
        self.pmSelector()
        self.parallaxSelector()
        self.importData()
        self.calculatingError()
        self.getCandidates()
        self.plottingCandidates()