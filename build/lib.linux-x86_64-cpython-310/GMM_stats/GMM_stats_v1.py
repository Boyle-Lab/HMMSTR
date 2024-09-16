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

from HMMSTR_utils.HMMSTR_utils import remove_outlier_IQR, throw_low_cov

class GMMStats:
    def __init__(self, target_row):
        self.name = target_row[0]
        self.prefix = target_row[1].rstrip().upper()[-30:]
        self.repeat = target_row[2].rstrip().upper()
        self.suffix = target_row[3].rstrip().upper()[:30]
    
    def plot_stats(self,counts_data,out):
        '''
        plot the histograms and median likelihood ratio line plots (hopefully with error bars)

        '''
        #(1) unmodified histogram
        to_plot = counts_data[['counts','freq']].drop_duplicates().sort_values(by="counts")
        if to_plot.shape[0] == 0 or to_plot.shape[1] == 0:
            return
        #(1) unmodified histogram
        sns.reset_orig()
        to_plot.plot.bar(x="counts",y="freq",figsize=(15,5),title= str(self.name), color="red", fontsize=7)
        plt.xlabel("Repeat count")
        plt.ylabel("Number of reads")
        fig = plt.gcf()
        fig.savefig(str(out +self.name+ "_supporting_reads_hist.pdf"))
        plt.show()
        plt.clf()
        return #currently does not work and isnt called
    def get_stats(self, out_count_file, out, plot_hists, filter_outliers=False):
        '''
        Function to calculate summary stats for a single target

        Args:
            out_count_name(str): output suffix for final result tsv

        Returns:
            counts_data(pandas dataframe): final calculations for given count file

        '''
        #read output file
        counts_data = pd.read_csv(out_count_file, sep=" ",header=None)
        counts_data.columns = ["read_id","strand","align_score","neg_log_likelihood","subset_likelihood","repeat_likelihood","align_start", "align_end","counts"]
        counts_data['freq'] = counts_data.groupby(by='counts')['read_id'].transform('count')
        #ADDED outlier filtering
        outliers = []
        if filter_outliers:
            filtered, outliers = remove_outlier_IQR(counts_data.counts)
            #counts_data = counts_data[counts_data.counts.isin(filtered)]
        counts_data['outlier'] = np.where(counts_data.counts.isin(outliers),True , False)
        #testing out including plotting here, will need to add a flag
        if plot_hists:
            self.plot_stats(counts_data, out)
        return counts_data

    def call_peaks(self,X_orig,out, num_max_peaks, plot=True, save_allele_plots = True):
        # https://www.cbrinton.net/ECE20875-2020-Spring/W11/gmms_notebook.pdf
        rng = np.random.RandomState(seed=1)

        #get transformed counts column to represent weighted frequencies
        X = X_orig.copy()
        X = X[["counts"]]
        if X.shape[0] < 2:
            curr_dict = {"name":self.name}
            #revised
            #check if we have 1 read
            if X.shape[0] == 1:
                # curr_dict["H1:mean"] = X.counts[0]
                curr_dict["H1:median"] = X.counts[0]
                curr_dict["H1:mode"] = X.counts[0]
                curr_dict["H1:SD"] = 0
                curr_dict["H1:supporting_reads"] = 1
                cluster_assignments = np.array([0])
            else:
                cluster_assignments = np.array([-1])
            for i in range(1,num_max_peaks+1):
                #curr_mean = "H"+str(i)+":mean"
                curr_median = "H"+str(i)+":median"
                curr_mode = "H"+str(i)+":mode"
                curr_sd = "H" + str(i) + ":SD"
                curr_support = "H" + str(i) + ":supporting_reads"
                #if curr_mean not in curr_dict.keys():
                   # curr_dict[curr_mean] = 0 #don't include the mean since our distributions are left skewed, report median and mode
                if curr_median not in curr_dict.keys():
                    curr_dict[curr_median] = 0
                if curr_mode not in curr_dict.keys():
                    curr_dict[curr_mode] = 0
                if curr_sd not in curr_dict.keys():
                    curr_dict[curr_sd] = 0
                if curr_support not in curr_dict.keys():
                    curr_dict[curr_support] = 0
        # for i in range(1,num_max_peaks+1):
        #     curr_dict["H"+str(i)+":SD"] = final_peaks_SD[i-1]
            curr_dict["num_supporting_reads"] = len(X_orig.index)
            return pd.Series(curr_dict), cluster_assignments
        # Fit models with 1-10 components
        #account for more components than samples
        if X.shape[0] < num_max_peaks: # FIXME this may be an issue in consistency, maybe save original number for output
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
                sns.set_style('whitegrid',{'grid.color': '#dcddde'}) # nice and clean grid
                if ax is None:
                    ax = plt.gca()

            # Compute PDF of whole mixture
            x = np.linspace(X.min(), X.max(), 1000)
            logprob = gmm.score_samples(x.reshape(-1, 1))
            pdf = np.exp(logprob)

            # Compute PDF for each component
            responsibilities = gmm.predict_proba(x.reshape(-1, 1))

            pdf_individual = responsibilities * pdf[:, np.newaxis]
            #print(pdf_individual.shape)
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
                    #ax.axvline(x = X.counts.mean(), linewidth=2, color='g', linestyle="dashed", label = "Mean")
                    ax.axvline(x = X.counts.median(), linewidth=2, color='purple',linestyle="dashed",label="Median")
                    #print("mean", X.counts.mean(), " mode: ", st.mode(np.array(X.counts))[0][0], "median: ", X.counts.median())
                if show_legend:
                    ax.legend()
                if not subplot: #save full peaks plot
                    plt.savefig(out + self.name + "peaks.pdf",bbox_inches='tight')
                    plt.show()
                    plt.clf()
                else:
                    plt.savefig(out + self.name + "allele_"+ str(i+1)+".pdf",bbox_inches='tight')
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
            plt.savefig(out + self.name + "AIC_BIC.pdf",bbox_inches='tight')
            plt.show()
            plt.clf()
            fig = plt.gcf()
            fig.set_size_inches(15,7)
        gmm_best = models[np.argmin(AIC)]
        max_peaks, ranks_idx = plot_mixture(gmm_best, X, title=self.name +": Gaussian Mixture Model of Best Fit")
        #print(max_peaks)
        # if plot:
        #     plt.title(self.name)
        #     plt.savefig(out + self.name + "peaks.pdf",bbox_inches='tight')
        #     plt.show()
        #     plt.clf()

        #ADDED 9/7/22 TO TEST CLUSTERING + NEW FITTING WITH SCIPY
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
                #plt.title(self.name)
                #plt.savefig(out + self.name + "allele_"+ str(i+1)+".pdf",bbox_inches='tight')
                #plt.show()
                #plt.clf()
            #add metrics to dictionary so they are outputted
            #cluster_stats[i]['mean'] = cluster.counts.mean()
            cluster_stats[i]['median'] = cluster.counts.median()
            cluster_stats[i]['mode'] = st.mode(np.array(cluster),keepdims=True)[0][0][0]
            if math.isnan(cluster.counts.std()):
                cluster_stats[i]['SD']  = 0
            else:
                cluster_stats[i]['SD'] = cluster.counts.std()
            cluster_stats[i]['supporting_reads'] = len(np.where(cluster_assignments == i)[0])
        #unravel cluster_stats dictionary for curr_row
        for i in cluster_stats.keys():
            #curr_dict["H" + str(i+1) + ":mean"] = cluster_stats[i]["mean"]
            curr_dict["H" + str(i+1) + ":median"] = cluster_stats[i]["median"]
            curr_dict["H" + str(i+1) + ":mode"] = cluster_stats[i]["mode"]
            curr_dict["H" + str(i+1) + ":SD"] = cluster_stats[i]["SD"]
            curr_dict["H" + str(i+1) + ":supporting_reads"] = cluster_stats[i]["supporting_reads"]

        #check for missing columns and fill them if not there
        for i in range(1,num_max_peaks+1):
            #curr_mean = "H"+str(i)+":mean"
            curr_median = "H"+str(i)+":median"
            curr_mode = "H"+str(i)+":mode"
            curr_sd = "H" + str(i) + ":SD"
            curr_support = "H" + str(i) + ":supporting_reads"
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
        # print(curr_row, "length: ", len(curr_row))
        # print(cluster_assignments)
        return curr_row, cluster_assignments
    
    def bootstrap_gmm(self,curr_row,data, resample_size, CI_width, max_peaks, out):
        #FIXME need to determine if I should only construct CIs from each allele separately instead of like how I am since if there is a large imbalance in supporting reads, the CIs will include both alleles by chance
        #throw out 1 read coverage, this is something to assume will happen if bootstrapping
        filtered = throw_low_cov(data)
        call_lists = {}
        for i in range(1,max_peaks+1):
            #call_lists["H"+str(i)+"_mean"]=[]
            call_lists["H"+str(i)+"_median"]=[]
            #call_lists["H"+str(i)+"_mode"]=[]

        for i in range(resample_size):
            #sample up to the depth of the experiment with replacement
            curr_sample = filtered.sample(frac=1,axis=0, replace = True)
            #call gmm with original freqs on this subset
            curr_row = self.call_peaks(curr_sample, out ,max_peaks, plot=False) 
            #add current call to lists
            for j in range(1,max_peaks+1):
                #call_lists["H"+str(j)+"_mean"].append(curr_row[0]['H'+str(j)+":mean"])
                call_lists["H"+str(j)+"_median"].append(curr_row[0]['H'+str(j)+":median"])
                #call_lists["H"+str(j)+"_mode"].append(curr_row[0]['H'+str(j)+":mode"])
        #mean_CIs = {}
        median_CIs = {}
        #mode_CIs = {}
        adjustment = float((1.0 - CI_width)/2)
        for j in range(1,max_peaks+1): #get confidence interval for each allele called
            #curr_means = call_lists["H"+str(j)+"_mean"]
            curr_medians = call_lists["H"+str(j)+"_median"]
            #curr_modes = call_lists["H"+str(j)+"_mode"]
            #curr_row[0]["H"+str(j)+":mean_CI"] = (np.quantile(curr_means,adjustment),np.quantile(curr_means,1-adjustment))
            curr_row[0]["H"+str(j)+":median_CI"] = (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            #curr_row[0]["H"+str(j)+":mode_CI"] = (np.quantile(curr_modes,adjustment),np.quantile(curr_modes,1-adjustment))
        #print(curr_row)
        return curr_row

    def bootstrap_gmm_allele_specific(self,data, resample_size, CI_width, out):
            #FIXME need to determine if I should only construct CIs from each allele separately instead of like how I am since if there is a large imbalance in supporting reads, the CIs will include both alleles by chance
            #throw out 1 read coverage, this is something to assume will happen if bootstrapping
            filtered = throw_low_cov(data)
            call_lists = {}
            #call_lists["mean"]=[]
            call_lists["median"]=[]
            #call_lists["mode"]=[]
    
            for i in range(resample_size):
                #sample up to the depth of the experiment with replacement
                curr_sample = filtered.sample(frac=1,axis=0, replace = True)
                #call gmm with original freqs on this subset
                curr_row = self.call_peaks(curr_sample, out ,1, plot=False) #single allele bootstrapping
                #add current call to lists
                #call_lists["mean"].append(curr_row[0]["H1:mean"])
                call_lists["median"].append(curr_row[0]["H1:median"])
                #call_lists["mode"].append(curr_row[0]["H1:mode"])
            #mean_CIs = {}
            median_CIs = {}
            #mode_CIs = {}
            adjustment = float((1.0 - CI_width)/2)
            #curr_means = call_lists["mean"]
            curr_medians = call_lists["median"]
            #curr_modes = call_lists["mode"]
            #mean_CI = (np.quantile(curr_means,adjustment),np.quantile(curr_means,1-adjustment))
            median_CI= (np.quantile(curr_medians,adjustment),np.quantile(curr_medians,1-adjustment))
            #mode_CI = (np.quantile(curr_modes,adjustment),np.quantile(curr_modes,1-adjustment))
            return median_CI