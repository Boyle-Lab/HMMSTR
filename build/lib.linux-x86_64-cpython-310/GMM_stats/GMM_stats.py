import seaborn as sns
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import scipy.stats as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
import math
from matplotlib.ticker import MaxNLocator
import warnings
warnings.filterwarnings('ignore')

from HMMSTR_utils.HMMSTR_utils import remove_outlier_IQR, throw_low_cov, remove_flanking_outliers

class GMMStats:
    def __init__(self, target_row):
        self.name = target_row[0]
        self.prefix = target_row[1].rstrip().upper()[-30:]
        self.repeat = target_row[2].rstrip().upper()
        self.suffix = target_row[3].rstrip().upper()[:30]
    
    def plot_stats(self,counts_data,out, strand=None):
        '''
        plot the read support histograms

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
        if strand is not None:
            fig.savefig(str(out + "_plots/"+self.name+ "_"+strand+"_supporting_reads_hist.pdf"))
        else:
            fig.savefig(str(out + "_plots/"+self.name+ "_supporting_reads_hist.pdf"))
        plt.show()
        plt.clf()
        return #currently does not work and isnt called
    def get_stats(self, out_count_file, out, plot_hists, filter_outliers=False, filter_quantile=0.25, flanking_like_filter=False, curr_strand=None):
        '''
        Function to calculate summary stats for a single target

        Args:
        -----------------------------------------------------------------------------
        out_count_name(str): output suffix for final result tsv
        out (str): output prefix
        plot_hists (bool): boolean indicating if histograms should be saved
        filter_outliers (bool): boolean indicating if outliers should be included or not
        filter_quantile (float): quantile that will be filtered if filter_outliers is true
        flanking_like_filter (bool): boolean designating if reads should be filtered based on the likelihood of their flanking sequence
        curr_strand (str): which strand we are performing operations on if applicable

        Returns:
        -----------------------------------------------------------------------------
        counts_data(pandas dataframe): final calculations for given count file

        '''
        #read output file
        counts_data = pd.read_csv(out_count_file, sep=" ",header=None)
        counts_data.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","repeat_start","repeat_end","align_start", "align_end","counts"]
        
        if curr_strand is not None:
            counts_data = counts_data[counts_data.strand == curr_strand]
        
        counts_data['freq'] = counts_data.groupby(by='counts')['read_id'].transform('count')
        #outlier filtering
        outliers = []
        flanking_outliers = []
        if filter_outliers:
            filtered, outliers = remove_outlier_IQR(counts_data.counts, filter_quantile)
        counts_data['outlier'] = np.where(counts_data.counts.isin(outliers),True , False)
        if flanking_like_filter:
            flanking_filtered, flanking_outliers = remove_flanking_outliers(counts_data)
        counts_data['flanking_outlier'] = np.where(counts_data.read_id.isin(flanking_outliers),True , False)
        if plot_hists:
            self.plot_stats(counts_data, out, curr_strand)
        return counts_data

    def call_peaks(self,X_orig,out, num_max_peaks, plot=True, save_allele_plots = True, strand=None):
        '''
        Call peaks using a Gaussian mixture model. Code adapted from https://www.cbrinton.net/ECE20875-2020-Spring/W11/gmms_notebook.pdf

        Args:
        -----------------------------------------------------------------------------
        X_orig (dataframe): counts matrix
        out (str): output prefix
        num_max_peaks (int): maximum number of peaks to call
        plot (bool): boolean indicating if plots should be saved
        save_allele_plots (bool): boolean indicating if allele specific plots should be saved
        strand (str): which strand operations are performed on if applicable

        Returns:
        -----------------------------------------------------------------------------
        curr_row (Series): series containing all stats for current target
        cluster_assignments (Series): series containing read to cluster assignments
        '''
        rng = np.random.RandomState(seed=1)

        #get transformed counts column to represent weighted frequencies
        X = X_orig.copy()
        X = X[["counts"]]
        if X.shape[0] < 2:
            curr_dict = {"name":self.name}
            #check if we have 1 read
            if X.shape[0] == 1:
                try:
                    # curr_dict["H1:mean"] = X.counts[0]
                    curr_dict["A1:median"] = X.counts[0]
                    curr_dict["A1:mode"] = X.counts[0]
                    curr_dict["A1:SD"] = 0
                    curr_dict["A1:supporting_reads"] = 1
                    cluster_assignments = np.array([0])
                except:
                    print("X.counts[0] does not exist, continuing to write null row...")
                    cluster_assignments = np.array([-1])
            else:
                cluster_assignments = np.array([-1])
            for i in range(1,num_max_peaks+1):
                #curr_mean = "A"+str(i)+":mean"
                curr_median = "A"+str(i)+":median"
                curr_mode = "A"+str(i)+":mode"
                curr_sd = "A" + str(i) + ":SD"
                curr_support = "A" + str(i) + ":supporting_reads"
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
            curr_dict["num_supporting_reads"] = len(X_orig.index)
            return pd.Series(curr_dict), cluster_assignments
        #account for more components than samples
        if X.shape[0] < num_max_peaks:
            num_max_peaks = X.shape[0]

        k_arr = np.arange(num_max_peaks) + 1
        models = [
        GaussianMixture(n_components=k).fit(X)
        for k in k_arr
        ]
        dfs = [] #subsets of X corresponding to each peak

        # Plot function
        def plot_mixture(gmm, X, show_legend=True, ax=None, title=None, subplot=False,cluster=None):
            if plot:
                plt.rcParams["figure.figsize"] = (10,6)
                sns.set_style('whitegrid',{'grid.color': '#dcddde'}) 
                if ax is None:
                    ax = plt.gca()

            # Compute PDF of whole mixture
            x = np.linspace(X.min(), X.max(), 1000)
            logprob = gmm.score_samples(x.reshape(-1, 1))
            pdf = np.exp(logprob)

            # Compute PDF for each component
            responsibilities = gmm.predict_proba(x.reshape(-1, 1))

            pdf_individual = responsibilities * pdf[:, np.newaxis]

            #get the max of each pdf
            maxInColumns = np.amax(pdf_individual, axis=0)

            # Plot data histogram
            if plot:
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                hist = ax.hist(X.counts, density=True, histtype='bar',  label='Data', bins=range(int(min(X.counts)),int(max(X.counts))+2),facecolor = '#9bc1cc', edgecolor='whitesmoke', linewidth=0.5, align="left") #updated plotting parameters
                ax.set_xlabel("Repeat Copy Number")
                ax.set_ylabel("Density")
                # Plot PDF of whole model
                if not subplot:
                    ax.plot(x, pdf, '-k', label='Mixture PDF',linewidth=2)

                # Plot PDF of each component, divided so we can call this multiple times, once for full mixture and once for each cluster
                if subplot:
                    ax.plot(x, pdf_individual[:,cluster],color="black")
                else:
                    ax.plot(x, pdf_individual, '--', label='Component PDF')
                if title is not None:
                    ax.set_title(title)
                if subplot:
                    #add metrics to plot
                    ymin, ymax = ax.get_ylim()
                    ax.set_ylim(ymin, ymax)
                    ax.vlines(x = st.mode(np.array(X.counts), keepdims=True)[0][0],ymin=ymin, ymax=ymax,linewidth=2, color='r',linestyle="dashed", label="Mode")
                    ax.axvline(x = X.counts.median(), linewidth=2, color='purple',linestyle="dashed",label="Median")
                if show_legend:
                    ax.legend()
                if not subplot: #save full peaks plot
                    if strand is not None:
                        plt.savefig(out + "_plots/"+ self.name +"_"+strand+ "peaks.pdf",bbox_inches='tight')
                    else:
                        plt.savefig(out + "_plots/"+ self.name + "peaks.pdf",bbox_inches='tight')
                    plt.show()
                    plt.clf()
                else:
                    if strand is not None:
                        plt.savefig(out + "_plots/"+ self.name + "_"+ strand+"allele_"+ str(i+1)+".pdf",bbox_inches='tight')
                    else:
                        plt.savefig(out + "_plots/"+ self.name + "allele_"+ str(i+1)+".pdf",bbox_inches='tight')
                    plt.show()
                    plt.clf()
                plt.show()
                plt.clf()
            return gmm.means_[np.argsort(maxInColumns)],np.argsort(maxInColumns)
        AIC = [m.aic(X) for m in models]
        BIC = [m.bic(X) for m in models]

        # Plot these metrics
        if plot:
            fig = plt.gcf()
            fig.set_size_inches(15,7)
            plt.plot(k_arr, AIC, label='AIC')
            plt.plot(k_arr, BIC, label='BIC')
            plt.xlabel('Number of Components ($k$)')
            plt.legend()
            plt.tight_layout()
            plt.title(self.name+": AIC and BIC")
            if strand is not None:
                plt.savefig(out + "_plots/"+ self.name +"_"+strand+ "_AIC_BIC.pdf",bbox_inches='tight')
            else:
                plt.savefig(out + "_plots/"+ self.name +"_AIC_BIC.pdf",bbox_inches='tight')
            plt.show()
            plt.clf()
            fig = plt.gcf()
            fig.set_size_inches(15,7)
        gmm_best = models[np.argmin(AIC)]
        max_peaks, ranks_idx = plot_mixture(gmm_best, X, title=self.name +": Gaussian Mixture Model of Best Fit")

        cluster_assignments = np.argmax(gmm_best.predict_proba(X), axis=1)
        for cluster in np.unique(cluster_assignments):
            dfs.append(X[np.where(cluster_assignments == cluster, True, False)])
        #plot separate clusters
        curr_dict = {"name":self.name}
        cluster_stats = {}
        for i,cluster in enumerate(dfs):
            #plot distributions separately
            cluster_stats[i] = {}
            if save_allele_plots and plot:
                test1, test2 = plot_mixture(gmm_best, cluster,title=self.name+": Model of Best Fit for Allele " + str(i+1), subplot=True,cluster=i)
            #add metrics to dictionary so they are outputted
            cluster_stats[i]['median'] = cluster.counts.median()
            cluster_stats[i]['mode'] = st.mode(np.array(cluster),keepdims=True)[0][0][0]
            if math.isnan(cluster.counts.std()):
                cluster_stats[i]['SD']  = 0
            else:
                cluster_stats[i]['SD'] = cluster.counts.std()
            cluster_stats[i]['supporting_reads'] = len(np.where(cluster_assignments == i)[0])
        #unravel cluster_stats dictionary for curr_row
        for i in cluster_stats.keys():
            #curr_dict["A" + str(i+1) + ":mean"] = cluster_stats[i]["mean"]
            curr_dict["A" + str(i+1) + ":median"] = cluster_stats[i]["median"]
            curr_dict["A" + str(i+1) + ":mode"] = cluster_stats[i]["mode"]
            curr_dict["A" + str(i+1) + ":SD"] = cluster_stats[i]["SD"]
            curr_dict["A" + str(i+1) + ":supporting_reads"] = cluster_stats[i]["supporting_reads"]

        #check for missing columns and fill them if not there
        for i in range(1,num_max_peaks+1):
            #curr_mean = "A"+str(i)+":mean"
            curr_median = "A"+str(i)+":median"
            curr_mode = "A"+str(i)+":mode"
            curr_sd = "A" + str(i) + ":SD"
            curr_support = "A" + str(i) + ":supporting_reads"
            #if curr_mean not in curr_dict.keys():
                #curr_dict[curr_mean] = 0
            if curr_median not in curr_dict.keys():
                curr_dict[curr_median] = 0
            if curr_mode not in curr_dict.keys():
                curr_dict[curr_mode] = 0
            if curr_sd not in curr_dict.keys():
                curr_dict[curr_sd] = 0
            if curr_support not in curr_dict.keys():
                curr_dict[curr_support] = 0
        curr_dict["num_supporting_reads"] = len(X_orig.index)
            
        curr_row = pd.Series(curr_dict)
        return curr_row, cluster_assignments
    
    def bootstrap_gmm(self,curr_row,data, resample_size, CI_width, max_peaks, out):
        '''
        Performs bootstrapping and constructs confidence intervals using GMM peakcalling

        Args:
        --------------------------------------------------------------------------------
        curr_row (series): current row in the genotype dataframe to append CIs to
        data (dataframe): counts dataframe to sample from
        resample_size (int): number of times to resample during bootstraping
        CI_width (float): width of confidience interval to calculate
        max_peaks (int): maximum number of peaks to call in each bootstrap iteration
        out (str): output prefix

        Returns:
        --------------------------------------------------------------------------------
        curr_row (series): the original input row with bootstrap results included
        '''
        #throw out 1 read coverage to avoid skewed bootstrapping results
        filtered = throw_low_cov(data)
        call_lists = {}
        for i in range(1,max_peaks+1):
            call_lists["A"+str(i)+"_median"]=[]

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
            curr_row = self.call_peaks(curr_sample, out ,max_peaks, plot=False) 
            #add current call to lists
            for j in range(1,max_peaks+1):
                call_lists["A"+str(j)+"_median"].append(curr_row[0]['A'+str(j)+":median"])
        adjustment = float((1.0 - CI_width)/2)
        for j in range(1,max_peaks+1): #get confidence interval for each allele called
            curr_medians = call_lists["A"+str(j)+"_median"]
            curr_row[0]["A"+str(j)+":median_CI"] = (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
        return curr_row[0]

    def bootstrap_gmm_allele_specific(self,data, resample_size, CI_width, out):
            '''
            Performs bootstrapping and constructs confidence interval for single allele using GMM peakcalling

            Args:
            --------------------------------------------------------------------------------
            data (dataframe): counts dataframe to sample from
            resample_size (int): number of times to resample during bootstraping
            CI_width (float): width of confidience interval to calculate
            out (str): output prefix

            Returns:
            --------------------------------------------------------------------------------
            median_CI (tuple): tuple containing the median confidence interval call for the current allele 
            '''
            #throw out 1 read coverage, this is something to assume will happen if bootstrapping
            filtered = throw_low_cov(data)
            call_lists = {}
            #call_lists["mean"]=[]
            call_lists["median"]=[]
            #call_lists["mode"]=[]

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
                curr_row = self.call_peaks(curr_sample, out ,1, plot=False) #single allele bootstrapping
                #add current call to lists
                call_lists["median"].append(curr_row[0]["A1:median"])
            adjustment = float((1.0 - CI_width)/2)
            curr_medians = call_lists["median"]
            median_CI= (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            return median_CI