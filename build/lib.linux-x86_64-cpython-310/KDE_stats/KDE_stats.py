import seaborn as sns
import scipy.stats as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
from numpy import array, linspace
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
import os
import statistics as stat
from collections import Counter

from HMMSTR_utils.HMMSTR_utils import throw_low_cov, remove_outlier_IQR, remove_flanking_outliers
warnings.filterwarnings('ignore')

class KDE_cluster:
    def __init__(self, target, out, out_count_name, discard_outliers,flanking_like_filter,curr_strand=None):
        self.name = target[0]
        self.out = out
        self.discard_outliers = discard_outliers #this will be used when we goto cluster
        self.flanking_like_filter = flanking_like_filter
        out_count_file = out + "_" + self.name + out_count_name
        if os.path.exists(out_count_file):
            counts_data = pd.read_csv(out_count_file, sep=" ",header=None)
            counts_data.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","repeat_start","repeat_end","align_start", "align_end","counts"]
            self.data = counts_data[counts_data.counts != 0] #discard 0 count reads
            if curr_strand is not None:
                self.strand = curr_strand
                self.data = counts_data[counts_data.counts != 0][counts_data.strand == curr_strand] #if this is a stranded run, we only want to keep data from the current strand
            else:
                self.strand = None
            if self.data.empty:
                self.data = None
        else:
            self.data = None
    def get_stats(self, plot_hists, filter_outliers=False, filter_quantile=0.25,flanking_like_filter=False, strand=None):
        '''
        Function to calculate summary stats for a single target

        Args:
        -----------------------------------------------------------------------------------------------
            plot_hists (bool): boolean indicating if read support histograms should be saved
            filter_outliers (bool): boolean indicating if outliers should be included in allele calls
            filter_quantile (float): if filter_outliers is set, this quantile of repeat copy number will be filtered
            flanking_like_filter (bool): if set, reads with excessively low likelihood in their flanking regions will be filtered
            strand (str): which strand to operated on if applicable

        Returns:
        -----------------------------------------------------------------------------------------------
            counts_data(pandas dataframe): final calculations for given count file

        '''
        #read output file
        counts_data = self.data
        counts_data['freq'] = counts_data.groupby(by='counts')['read_id'].transform('count')
        #outlier filtering
        outliers = []
        flanking_outliers=[]
        if filter_outliers:
            filtered, outliers = remove_outlier_IQR(counts_data.counts, filter_quantile)
            #counts_data = counts_data[counts_data.counts.isin(filtered)]
        counts_data['outlier'] = np.where(counts_data.counts.isin(outliers),True , False)
        #check for flanking likelihood outliers
        if flanking_like_filter:
            flanking_filtered, flanking_outliers = remove_flanking_outliers(counts_data)
        counts_data['flanking_outlier'] = np.where(counts_data.read_id.isin(flanking_outliers),True , False)
        if plot_hists:
            self.plot_stats(counts_data, self.out, strand)
        return counts_data
    def plot_stats(self,counts_data,out, strand=None):
        '''
        plot the histograms

        Args:
        -----------------------------------------------------------------------------
        counts_data (dataframe): dataframe continaing copy numbers from Viterbi output
        out (str): output prefix
        strand (str): strand we are plotting if applicable

        Returns:
        ----------------------------------------------------------------------------
        None

        '''
        to_plot = counts_data[['counts','freq']].drop_duplicates().sort_values(by="counts")
        if to_plot.shape[0] == 0 or to_plot.shape[1] == 0:
            return
        sns.reset_orig()
        to_plot.plot.bar(x="counts",y="freq",figsize=(15,5),title= str(self.name), color="red", fontsize=7)
        plt.xlabel("Repeat count")
        plt.ylabel("Number of reads")
        fig = plt.gcf()
        if self.strand is not None:
            fig.savefig(str(self.out + "_plots/"+self.name+ "_"+ self.strand +"_supporting_reads_hist.pdf"))
        else:
            fig.savefig(str(self.out + "_plots/"+self.name+"_supporting_reads_hist.pdf"))
        plt.show()
        plt.clf()
        return #currently does not work and isnt called
    def call_clusters(self,kernel="gaussian", bandwidth = 'scott', max_k = 2, output_plots = False, subset=None, filter_quantile=0.25, allele_specific = False, allele = None, strand=None):
        '''
        This function uses kernal density estimation to resolve alleles. It is ideal when data is normally distributed
        and has a relatively small range or overlapping count distributions

        Args
        --------------------------------------------------------------------------------------------------------------
        kernel: str
            String designating kernel to use for the kde
        bandwidth: float or str
            Bandwidth to use in kde
        max_k: int
            The maximum number of clusters to keep in outputs, all others are filtered by cluster size
        output_plots: bool
            True if density plot will be saved, False if not
        subset: bool
            if True, cluster a subset of the data
        filter_quantile: float
            qunatile to filter if filter outliers is set
        allele_specific: bool
            if True, plot allele specific plots
        allele: str
            allele to plot and cluster
        strand: str
            which strand to plot and cluster if applicable

        Returns
        -------------------------------------------------------------------------------------------------------------
        clusters: dict
            Dictionary containing all counts assigned to each cluster
        allele_calls: dict
            Dictionary containing the allele call from each cluster as well as the size and standard deviation within the cluster
        outliers: pandas Series
            Series of counts that are identified as outliers and were dicarded before clustering. Empty if filter_outliers is False
        flanking_outliers: pandas Series
            Series of reads to filter based on the likelihood of the viterbi path through the prefix and suffix

        '''
        #filter outliers if called for
        outliers = pd.Series()
        flanking_outliers = pd.Series()
        if self.discard_outliers and self.flanking_like_filter == False: #only discard outliers
            filtered, outliers = remove_outlier_IQR(self.data.counts,filter_quantile)
            a = np.array(filtered).reshape(-1,1)
        elif self.flanking_like_filter: #either flanking filter or both
            flanking_filtered, flanking_outliers = remove_flanking_outliers(self.data)
            #check if we also filtered outliers so we can do both
            if self.discard_outliers:
                filtered, outliers = remove_outlier_IQR(flanking_filtered.counts,filter_quantile) #filter outliers
                a = np.array(flanking_filtered[flanking_filtered.counts.isin(filtered.unique())].counts).reshape(-1,1) #use reads that are both in both of our filtered sets
            else:
                a = np.array(flanking_filtered.counts).reshape(-1,1)

        else:
            a = np.array(self.data.counts).reshape(-1,1)
        if subset is not None:
            a = np.array(subset.counts).reshape(-1,1) #allow a subset of data to be passed for reclustering
        try:#initial call
            kde = KernelDensity( bandwidth=bandwidth,kernel=kernel).fit(a)

            #trying minimum bandwidth
            if kde.bandwidth_ < 0.5:
                kde = KernelDensity( bandwidth=0.5,kernel=kernel).fit(a)
        except:
            if len(a) == 0:
                clusters = -1
                allele_calls = {"A"+str(k)+":median":-1 for k in range(1,max_k+1)}
                allele_calls.update({"A"+str(k)+":mode":-1 for k in range(1,max_k+1)})
                allele_calls.update({"A"+str(k)+":SD":-1 for k in range(1,max_k+1)})
                allele_calls.update({"A"+str(k)+":supporting_reads":-1 for k in range(1,max_k+1)})
                allele_calls["bandwidth"] = -1
                return clusters,allele_calls , outliers, flanking_outliers

            print("Switching to a constant bandwidth=0.5")
            kde = KernelDensity( bandwidth=0.5,kernel=kernel).fit(a)
        s = linspace(min(a)-2,max(a)+2)
        e = kde.score_samples(s.reshape(-1,1))
        mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

        #run and change bandwidth till # maxima matches or is less than max_k
        while len(mi)+1 > max_k:
            bandwidth = kde.bandwidth_ + 0.1
            kde = KernelDensity( bandwidth=bandwidth,kernel=kernel).fit(a)
            e = kde.score_samples(s.reshape(-1,1))
            mi, ma = argrelextrema(e, np.less)[0], argrelextrema(e, np.greater)[0]

        #save cluster info based on minima and maxima
        all_clusters = {}
        maxima = {}
        sizes = {}
        if len(mi) == 0:
            all_clusters["1"] = a
            maxima["1"] = s[ma]
            sizes["1"] = len(a)
        elif len(mi) == 1:
            #print("hi")
            all_clusters["1"] = a[a < s[mi][0]]
            maxima["1"] = s[ma][0]
            sizes["1"] = len(a[a < s[mi][0]])
            all_clusters["2"] = a[a >= s[mi][0]]
            if len(s[ma]) > 1:
                maxima["2"] = s[ma][1]
            else:
                if isinstance(stat.mode(all_clusters["2"]), np.float64):
                    maxima["2"] = stat.mode(all_clusters["2"])
                else:
                    maxima["2"] = stat.mode(all_clusters["2"])[0] #if multiple modes returned, pick the first
            sizes["2"] = len(a[a >= s[mi][0]])
        else: # more than 2 clusters
            all_clusters["1"] = a[a < s[mi][0]]
            maxima["1"] = s[ma][0]
            sizes["1"] = len(a[a < s[mi][0]])
            for i_cluster in range(len(mi)-1):
                all_clusters[str(i_cluster+2)] = a[(a >= s[mi][i_cluster]) * (a <= s[mi][i_cluster+1])]
                maxima[str(i_cluster+2)] = s[ma][i_cluster+1]
                sizes[str(i_cluster+2)] = len(a[(a >= s[mi][i_cluster]) * (a <= s[mi][i_cluster+1])])
            all_clusters[str(len(mi)+1)] = a[a >= s[mi][-1]]
            maxima[str(len(mi)+1)] = s[ma][-1]
            sizes[str(len(mi)+1)] = len(a[a >= s[mi][-1]])

        #filter clusters based on max_k
        sizes_series = pd.Series(sizes)
        #sort sizes
        sizes_series_sorted = sizes_series.sort_values(ascending=False)
        #get max_k clusters by size
        k_clusters = list(sizes_series_sorted.head(max_k).index)
        clusters = {k: np.unique(all_clusters[k]) for k in k_clusters}
        non_one_maxima = False
        for k in k_clusters:
            if isinstance(maxima[k], np.float64) == False and len(maxima[k]) != 1:
                non_one_maxima = True
        if non_one_maxima == False:
            allele_calls = {}
            allele_calls.update({"A" +str(k)+":supporting_reads":sizes[k] for k in k_clusters})
            allele_calls.update({"A" +str(k)+":SD":np.std(all_clusters[k]) for k in k_clusters})
            allele_calls.update({"A" +str(k)+":median":np.median(all_clusters[k]) for k in k_clusters})
            allele_calls.update({"A" +str(k)+":mode":int(st.mode(np.array(all_clusters[k]),keepdims=True)[0][0]) for k in k_clusters}) #I dont know if this will be the same, for some reason my results from margits data dont line up with finding the mode
            
            allele_calls["bandwidth"] = kde.bandwidth_
        elif len(k_clusters) == 1:
            #this is the case when there is only 1 value detected thus it is the only allele that can be called
            allele_calls = {}
            allele_calls.update({"A1:supporting_reads":sizes["1"]})
            allele_calls.update({"A1:SD":np.std(all_clusters["1"])})
            allele_calls.update({"A1:median":np.median(all_clusters["1"])})
            allele_calls.update({"A1:mode":int(st.mode(np.array(all_clusters["1"]),keepdims=True)[0][0])})
            allele_calls["bandwidth"] = kde.bandwidth_
        else:
            allele_calls = {}
            allele_calls["bandwidth"] = -1
            allele_calls.update({"A"+str(k)+":SD":0 for k in k_clusters})
            allele_calls.update({"A"+str(k)+":median":0 for k in k_clusters})
            allele_calls.update({"A"+str(k)+":mode":0 for k in k_clusters})
            return clusters,allele_calls , outliers, flanking_outliers

        if output_plots:
            if len(mi) == 0:
                try:
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                        
                    if len(s[ma] != 0):
                        plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                    
                except:
                    print("Encounted an exception in homozygous plotting case (KDE)")
                if allele_specific:
                    plt.title(self.name + " KDE: Allele " + allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name +"_" +strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_allele_"+ allele+"_KDE.pdf"))
                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name +"_" +strand+ "_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
            elif len(mi) == 1:
                try:
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                    plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                except:
                    print("An exception has occured")
                plt.plot(s[mi], np.exp(e[mi]), 'ro',label="KDE Minima")
                if allele_specific:
                    plt.title(self.name + " KDE: Allele "+ allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name +"_"+ strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_allele_"+ allele+"_KDE.pdf"))

                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name +"_"+ strand+ "_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
            #dynamically print clusters
            else:
                try:
                    plt.hist(a, density=True, histtype='bar',  label='Data', bins=range(int(min(a)),int(max(a))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left")
                    plt.plot(s[:],np.exp(e[:]),color="black",lw=2,linestyle="-",label="Kernel Density Estimate")
                    plt.plot(s[ma], np.exp(e[ma]), 'bo',label="KDE Maxima")
                except:
                    print("An exception occured in plotting the KDE")
                plt.plot(s[mi], np.exp(e[mi]), 'ro',label="KDE Minima")
                if allele_specific:
                    plt.title(self.name + " KDE: Allele " +allele)
                    ymin, ymax = plt.ylim()
                    plt.ylim(ymin, ymax)
                    plt.vlines(x = allele_calls["A1:mode"],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    plt.vlines(x = allele_calls["A1:median"], ymin=ymin, ymax=ymax,linewidth=2, color='purple',linestyle="dashed",label="Median")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_"+ strand+"_allele_"+ allele+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+self.name + "_allele_"+ allele+"_KDE.pdf"))

                else:
                    plt.title(self.name + " KDE")
                    plt.legend()  
                    fig = plt.gcf()
                    if strand is not None:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_"+ strand+"_KDE.pdf"))
                    else:
                        fig.savefig(str(self.out + "_plots/"+ self.name + "_KDE.pdf"))
                plt.show()
                plt.clf()
        return clusters, allele_calls, outliers, flanking_outliers
    def assign_clusters(self, clusters, outliers, flanking_outliers):
        '''
        Function to assign individual reads to clusters determined by KDE

        Parameters
        ------------------------------------------------------------------
        clusters: dict
            Dictionary of counts assigned to each cluster
        outliers: series
            series containing all outliers identified
        flanking_outliers: series
            series containing all flanking sequence likelihood outliers
        
        Returns:
        ----------------------------------------------------------------------
        assignemnt_df: pandas DataFrame
            DataFrame equivalent to input count data with additional columns for cluster
            assignemnts and if the count was an outlier
        '''
        count_assignemnts = {}
        if clusters == -1:
            return pd.DataFrame(columns =['read_id', 'strand', 'align_score', 'neg_log_likelihood', 'subset_likelihood', 'repeat_likelihood', 'repeat_start', 'repeat_end', 'align_start', 'align_end', 'counts', 'freq', 'cluster_assignments', 'outlier', 'flanking_outlier'])
        for cluster in clusters:
            for curr in clusters[cluster]:
                count_assignemnts[curr] = cluster
        labels = pd.DataFrame()
        labels['counts'] = count_assignemnts.keys()
        labels['cluster_assignments'] = count_assignemnts.values()
        assignment_df = self.data.merge(labels, how = 'left', on = 'counts')
        assignment_df['outlier'] = np.where(assignment_df.counts.isin(outliers),True , False)
        assignment_df['flanking_outlier'] = np.where(assignment_df.read_id.isin(flanking_outliers),True , False)
        return assignment_df

    def bootstrap_KDE(self,data, resample_size, CI_width, max_peaks, out):
        '''
        Performs bootstrapping and constructs confidence intervals using GMM peakcalling

        Args:
        --------------------------------------------------------------------------------
        data (dataframe): counts dataframe to sample from
        resample_size (int): number of times to resample during bootstraping
        CI_width (float): width of confidience interval to calculate
        max_peaks (int): maximum number of peaks to call in each bootstrap iteration
        out (str): output prefix

        Returns:
        --------------------------------------------------------------------------------
        curr_row (series): the original input row with bootstrap results included
        '''
        filtered = throw_low_cov(data)
        call_lists = {}
        for i in range(1,max_peaks+1):
            call_lists["A"+str(i)+":median"]=[]

        for i in range(resample_size):
            if resample_size > 100 and (i+1)%10 == 0:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            elif resample_size < 100 and resample_size > 10 and (i+1)%5 == 0:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            elif resample_size < 10:
                print("Round " + str(i+1)+ "/" + str(resample_size))
            #sample up to the depth of the experiment with replacement
            curr_sample = filtered.sample(frac=1,axis=0, replace = True)
            #call gmm with original freqs on this subset
            clusters, allele_calls, outliers, flanking_outliers = self.call_clusters(kernel="gaussian", bandwidth = 'scott', max_k = max_peaks, output_plots = False, subset=curr_sample, filter_quantile=0.25) 
            #add current call to lists
            for j in range(1,max_peaks+1):
                if ("A"+str(j)+":median") not in allele_calls.keys():
                    call_lists["A"+str(j)+":median"].append(0)
                    allele_calls["A"+str(j)+":median"] = 0
                    continue
                call_lists["A"+str(j)+":median"].append(allele_calls["A"+str(j)+":median"])
        median_CIs = {}
        adjustment = float((1.0 - CI_width)/2)
        for j in range(1,max_peaks+1): #get confidence interval for each allele called
            curr_medians = call_lists["A"+str(j)+":median"]
            median_CIs["A"+str(j)+":median_CI"] = (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
        return median_CIs

    def bootstrap_KDE_allele_specific(self,data, resample_size, CI_width, out):
            #throw out 1 read coverage, this is something to assume will happen if bootstrapping
            filtered = throw_low_cov(data)
            call_lists = {}
            call_lists["median"]=[]
    
            for i in range(resample_size):
                if resample_size > 100 and (i+1)%10 == 0:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                elif resample_size < 100 and resample_size > 10 and (i+1)%5 == 0:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                elif resample_size < 10:
                    print("Round " + str(i+1)+ "/" + str(resample_size))
                #sample up to the depth of the experiment with replacement
                curr_sample = filtered.sample(frac=1,axis=0, replace = True)
                #call gmm with original freqs on this subset
                clusters, allele_calls, outliers, flanking_outliers = self.call_clusters(kernel="gaussian", bandwidth = 'scott', max_k = 1, output_plots = False, subset=curr_sample, filter_quantile=0.25)  #single allele bootstrapping
                #add current call to lists
                call_lists["median"].append(allele_calls["A1:median"])
            adjustment = float((1.0 - CI_width)/2)
            curr_medians = call_lists["median"]
            median_CI= (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            return median_CI