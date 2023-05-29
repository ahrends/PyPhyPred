#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue, May 30th 2023
Written during WIN FMRIB/OHBA Analysis groups' PyPhyPred hackathon
#add names   
A separate copy of the licence is also included within the repository.

This file is TEMPORARY, includes some misc functions that might be useful at different parts of the package
 
"""

###############################################################################

#import python packages that you need here; e.g. numpy scipy sklearn etc
import numpy as np
import scipy
import copy 
###############################################################################

def exclude_keyword(keyword,conf_names_all):
    names_keyword=[aa for aa in conf_names_all if keyword in aa]
    include_keyword_inds=np.asarray([conf_names_all.index(aa) for aa in names_keyword])
    exclude_keyword_inds=np.asarray([aa for aa in range(len(conf_names_all)) if aa not in include_keyword_inds])
    conf_names_exclude=[copy.deepcopy(conf_names_all[ii]) for ii in exclude_keyword_inds]
    return include_keyword_inds,exclude_keyword_inds,conf_names_exclude

#names_site3=[aa for aa in conf_names_all if 'Site_3' in aa]
#site3_inds=np.asarray([conf_names_all.index(aa) for aa in names_site3])
#site12_inds=np.asarray([aa for aa in range(len(conf_names_all)) if aa not in site3_inds])
#conf_names_site12=[conf_names_all[ii] for ii in site12_inds]
    
    
def partial_corr_precmat(C):
    covariance_estimator = skcov.GraphLasso(alpha=0,verbose=1) #choose 0<alpha<1 to get a sparse covariance/precision
    covariance_estimator.fit(C)
    precmat=covariance_estimator.get_precision()
    Cstd = np.sqrt(np.diag(precmat))
    #Cstd[Cstd < 1.0e-5] = 1.0
    return precmat / np.outer(Cstd, Cstd)

def normalise_precmat(precmat):
    Cstd = np.sqrt(np.diag(precmat))
    #    Cstd[Cstd < 1.0e-5] = 1.0
    return precmat / np.outer(Cstd, Cstd)


def demean_mat(X,dim):
    if np.min(X.shape)<=1:
        Xm=np.mean(X)
    else:
        if dim==0:
            Xm=np.tile(np.mean(X,dim)[np.newaxis,:],(X.shape[0],1))
        else:
            Xm=np.tile(np.mean(X,dim)[:,np.newaxis],(1,X.shape[1]))
    Xdem=X-Xm
    return Xdem

def normalise_by_var_mat(X,dim=0,global_norm=False):
    if global_norm:
        Xnorm=X/np.std(X)
    else:
        if dim==0:
            Xnorm=X/np.tile(np.std(X,dim)[np.newaxis,:],(X.shape[0],1))
        else:
            Xnorm=X/np.tile(np.std(X,dim)[:,np.newaxis],(1,X.shape[1]))
    return Xnorm

def binarise_mat(X,thresh=0.5,num_nonzero=None):
    if num_nonzero is None:
        X_bin=X.copy()
        X_bin[np.abs(X)<thresh]=0
        X_bin[np.abs(X_bin)>0]=1
    else:
        X_bin=np.zeros(X.shape)
        for coli,col_num in enumerate(num_nonzero):
            thiscol=X[:,coli].copy()
            thislocs=np.argsort(np.abs(thiscol))[::-1][:col_num]
            X_bin[thislocs,coli]=1
    return X_bin


def gaussianise_nonparametric(Xmat,nquantiles=None,X_trans_qt=None):
    #Xmat should be nsample x nfeature: it will be gaussianised "per column"
    Xmat_tmp=Xmat.copy()
    if nquantiles is None:
        nquantiles=int(Xmat_tmp.shape[0]/5)
    if X_trans_qt is None:
        qt = QuantileTransformer(n_quantiles=nquantiles, output_distribution='normal',random_state=0)    
        X_trans_qt = qt.fit(Xmat_tmp)
    Xmat_normal=X_trans_qt.transform(Xmat_tmp)
    Xmat_normal = np.nan_to_num(Xmat_normal)
    return Xmat_normal,X_trans_qt

def compute_svd_basic(raw_mat,svd_perc=0.9,N_svd=None,demean_cols=True):
    if demean_cols:
        raw_mat1=demean_mat(raw_mat,0)
    else:
        raw_mat1=raw_mat.copy()
        #raw_mat1 = np.nan_to_num(raw_mat1)
    
    U_svd,s_svd,V_svd = np.linalg.svd(raw_mat1, full_matrices=False)
    if N_svd is None:
        N_svd=np.where(np.cumsum(s_svd**2/np.sum(s_svd**2))>svd_perc)[0][0]
    else:
        N_svd_cal=np.where(np.cumsum(s_svd**2/np.sum(s_svd**2))>svd_perc)[0][0]
        print(str(svd_perc)+" of variance is explained by "+str(N_svd_cal)+" PCs")
    svd_transformed_mat = U_svd[:,:N_svd] @ np.diag(s_svd[:N_svd]) 
    return svd_transformed_mat

def deconfound_mat(Xmat,confmat,train_beta=None):
    Xmat_deconf=demean_mat(Xmat,0)
    conf_mat=demean_mat(confmat,0)
    if train_beta is None:
        train_beta=scipy.linalg.pinv(conf_mat)@Xmat_deconf
    Xmat_deconf=Xmat_deconf-conf_mat@train_beta
    Xmat_deconf=demean_mat(Xmat_deconf,0)
    return Xmat_deconf ,train_beta

def deconfound_mat_fslnets(Xmat,confmat,train_beta=None):
    Xmat_deconf=demean_mat(Xmat,0)
    conf_mat=demean_mat(confmat,0)   
    Xmatd=Xmat_deconf.copy()#Xmats is yd; Xmat_deconf is y
    if not np.isnan(conf_mat).any():
        #conf_mat=compute_svd_basic(conf_mat,N_svd=np.linalg.matrix_rank(conf_mat),demean_cols=False)
        #print(conf_mat.shape)
        Xmatd[np.isnan(Xmatd)]=0 #this is allowed because setting NaNs to zero (in Xmat) doesn't change the maths...much....
        if train_beta is None:
            #print(scipy.linalg.pinv(conf_mat).shape)
            train_beta=scipy.linalg.pinv(conf_mat)@Xmatd
            train_beta[np.abs(train_beta)<1e-10]=0
        Xmatd=Xmat_deconf-conf_mat@train_beta # re-inherit nans
        Xmatd=demean_mat(Xmatd,0)
    else:        
        Xmatd[:]=np.nan#set to all NaN because we are not going to necessarily write into all elements below
        for cnt in range(Xmatd.shape[1]):
            no_nan_inds=np.isnan(np.sum(conf_mat,1)+Xmat_deconf[:,cnt])==False
            conf_mat_nonan=conf_mat[no_nan_inds,:].copy()
            conf_mat_nonan=demean_mat(conf_mat_nonan,0)
            #conf_mat_nonan=compute_svd_basic(conf_mat_nonan,N_svd=np.linalg.matrix_rank(conf_mat_nonan),demean_cols=False)
            if train_beta is None:
                train_beta=scipy.linalg.pinv(conf_mat_nonan)@Xmat_deconf[no_nan_inds,cnt]
                train_beta[np.abs(train_beta)<1e-10]=0
            Xmatd[no_nan_inds,cnt]=Xmat_deconf[no_nan_inds,cnt]-conf_mat_nonan@train_beta
            Xmatd[no_nan_inds,cnt]=demean_mat(Xmatd[no_nan_inds,cnt],0)
        
    return Xmatd ,train_beta 


def cross_validated_regression(Xmat_train,Xmat_test,target_train,target_type='continuous',options={}):
    options_local=copy.deepcopy(options)
    start_time = time.time()
    if target_type=='binary' or target_type=='multinomial':#binary or multinomial- do logistic regression
        print("binary or multinomial target")
        # all_keys={'Cs':5, 'cv':5, 'penalty':'elasticnet', 'solver':'saga', 'tol':0.01, 'max_iter':100, 'l1_ratios':[.1, .5, .7, .9, .95, .99, 1],'multi_class':'ovr','n_jobs':3}
        # for this_key in list(all_keys.keys()):
        #     if this_key not in options_local.keys():
        #         options_local[this_key] = copy.deepcopy(all_keys[this_key])
        
        # regr=LogisticRegressionCV(Cs=options_local['Cs'], cv=options_local['cv'], penalty=options_local['penalty'], solver=options_local['solver'], tol=options_local['tol'], l1_ratios=options_local['l1_ratios'],multi_class=options_local['multi_class'],n_jobs=options_local['n_jobs'])
        all_keys={'C':np.logspace(-4,4,5), 'cv':5, 'penalty':'elasticnet', 'solver':'saga', 'tol':0.01, 'max_iter':100, 'l1_ratio':[.1, .5, .7, .9, .95, .99, 1],'multi_class':'ovr','n_jobs':3,'n_iter':10}
        for this_key in list(all_keys.keys()):
            if this_key not in options_local.keys():
                options_local[this_key] = copy.deepcopy(all_keys[this_key])
        print(options_local)
        start_time = time.time()
        regr_pre = LogisticRegression(penalty=options_local['penalty'], solver=options_local['solver'], tol=options_local['tol'],multi_class=options_local['multi_class'])
        parameters = {'l1_ratio':options_local['l1_ratio'], 'C':options_local['C']}
        regr = RandomizedSearchCV(regr_pre, parameters,cv=options_local['cv'],n_jobs=options_local['n_jobs'],n_iter=options_local['n_iter'])

    elif target_type=='continuous':#continuous- do linear cross-validated elasticnet
        print("continuous target")
        all_keys={'l1_ratio':[0.1,0.5,0.7,0.9,0.95,0.99,1.0], 'cv':5, 'n_alphas':10, 'kvnormalize':False,'eps':1/(10.0*10),'tol':0.001,'n_jobs':None}
        for this_key in list(all_keys.keys()):
            if this_key not in options_local.keys():
                options_local[this_key] = copy.deepcopy(all_keys[this_key])
        regr=ElasticNetCV(l1_ratio=options_local['l1_ratio'], n_alphas=options_local['n_alphas'], cv=options_local['cv'], normalize=options_local['kvnormalize'], tol=options_local['tol'], eps=options_local['eps'],n_jobs=options_local['n_jobs'])
    elif target_type=='count':#positive counts- do poisson regression
        print("count integer target")
        all_keys={'reg_lambda':0.8*np.exp(-0.4 * np.arange(25)), 'distr':['poisson'], 'alpha':[.1, .5, .7, .9, .95, .99, 1],'tol':[0.01],'cv':5,'n_jobs':5,'n_iter':10,'max_iter': [100]}
        for this_key in list(all_keys.keys()):
            if this_key not in options_local.keys():
                options_local[this_key] = copy.deepcopy(all_keys[this_key])
        parameters = {'alpha': options_local['alpha'],
                       'reg_lambda': options_local['reg_lambda'],
                       'distr':options_local['distr'],#'gaussian',
                       'max_iter': options_local['max_iter'],
                       'tol': options_local['tol']}
        regr_pre=GLM()
        # regr = RandomizedSearchCV(MyGlmnetCV(), parameters,cv=options_local['cv'],n_jobs=options_local['n_jobs'],n_iter=options_local['n_iter'])
        regr = GridSearchCV(regr_pre, parameters,cv=options_local['cv'],n_jobs=options_local['n_jobs'])#,n_iter=options_local['n_iter'])

        # all_keys={'lambdau':np.asarray([[thislam] for thislam in 0.8*np.exp(-0.2 * np.arange(50))]), 'cv':5,  'family':['poisson'], 'alpha':[.1, .5, .7, .9, .95, .99, 1],'n_jobs':5,'n_iter':10,'max_iter': [100]}
        # for this_key in list(all_keys.keys()):
        #     if this_key not in options_local.keys():
        #         options_local[this_key] = copy.deepcopy(all_keys[this_key])
        # parameters = {'alpha': options_local['alpha'],
        #                'lambdau': options_local['lambdau'],
        #                'family':options_local['family'],#'gaussian',
        #                'max_iter': options_local['max_iter']}
        # # regr = RandomizedSearchCV(MyGlmnetCV(), parameters,cv=options_local['cv'],n_jobs=options_local['n_jobs'],n_iter=options_local['n_iter'])
        # regr = GridSearchCV(MyGlmnetCV(), parameters,cv=options_local['cv'],n_jobs=options_local['n_jobs'])#,n_iter=options_local['n_iter'])
    elif target_type=='ordinal':
        print("ordinal target, doing regression using ElasticNetCV for now but this has to change")
        all_keys={'l1_ratio':[0.1,0.5,0.7,0.9,0.95,0.99,1.0], 'cv':5, 'n_alphas':10, 'kvnormalize':False,'eps':1/(10.0*10),'tol':0.001,'n_jobs':None}
        for this_key in list(all_keys.keys()):
            if this_key not in options_local.keys():
                options_local[this_key] = copy.deepcopy(all_keys[this_key])
        regr=ElasticNetCV(l1_ratio=options_local['l1_ratio'], n_alphas=options_local['n_alphas'], cv=options_local['cv'], normalize=options_local['kvnormalize'], tol=options_local['tol'], eps=options_local['eps'],n_jobs=options_local['n_jobs'])
    else:
        raise ValueError('target type unknown: '+target_type)
    
    regr.fit(demean_mat(Xmat_train,0), target_train)
    predtarget=regr.predict(demean_mat(Xmat_test,0))
    print("--- %s seconds ---" % (time.time() - start_time))
    return predtarget 

def density_scatter(x,y,plot_baseline=True,plot_xy=False,fig=None,ax=None):
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    ax.scatter(x, y, c=z, s=50)#, edgecolor='')
    if plot_baseline:
        ax.plot(x,np.zeros(x.shape),'--r')
    if plot_xy:
        ax.plot(x,x,'--r')
    #plt.show(block=False)
    return ax

def cross_validation_main(Xmat_pre,conf_mat_raw,sorted_subjects,thisclass_path_stoch,dict_target_vectors,dict_target_types,
                          global_kvnums,required_target,svd_feature_selection=False,svd_perc_explainedvar=0.95,xvalidated_featselect=False,
                          xvalidated_featselect_perc=None,xvalidated_featselect_num=None,gaussianise_dists=[True,True,True],kvnum=5,global_kvnum=None,global_deconf=False,
                          target_deconf=False,xmat_deconf=True,nodeconf_predict=False,cv_options={},do_subjICA=False,subjICA_prefix='',required_IDP='',minimal_conf=False,name_prefix_conf="_minConf"):

    """
    some x-validation arguments
    """  
    if global_kvnum is None:
        global_kvnum=len(global_kvnums)    
    
    
    name_prefix=""
    if svd_feature_selection:
        name_prefix=name_prefix+"SVD_pnt"+str(svd_perc_explainedvar)[2:]
    if xvalidated_featselect:
        if xvalidated_featselect_perc is not None:
            name_prefix=name_prefix+"corr_pnt"+str(xvalidated_featselect_perc)[2:]
        else:
            name_prefix=name_prefix+"corr_num"+str(xvalidated_featselect_num)
    if not svd_feature_selection and not xvalidated_featselect:
        name_prefix=name_prefix+"raw"
    
    if minimal_conf:
        name_prefix=name_prefix+name_prefix_conf#"_minConf"

    
    
    if do_subjICA:
        name_prefix=name_prefix+subjICA_prefix#"_subjICA"#+required_IDP
    print(name_prefix)
    
    """
    if you need dimensionality reduction outside x-validation loop
    """
    
    if svd_feature_selection:
        Xmat_pre=compute_svd_basic(Xmat_pre,svd_perc=svd_perc_explainedvar)
    else:
        N_svd=Xmat_pre.shape[1]#
    
    #selected_features =np.arange(N_svd)  
    """
    create feature matrix, globally deconfounded feature matrix and target variables (i.e. age and sex)
    """
    #SRW_stoch=demean_mat(SRW_stoch,0)
    ##sgcorrstoch_raw=SRW_stoch.copy()
    ##SRW_stoch=SRW_stoch-np.tile(np.mean(SRW_stoch,1)[:,np.newaxis],(1,SRW_stoch.shape[1]))
    #SRW_stoch=SRW_stoch/np.max(np.abs(SRW_stoch))
    #
    ##SRW_stoch = scipy.stats.mstats.zscore(SRW_stoch, axis=1)
    #
    #SRW_stoch = np.nan_to_num(SRW_stoch)
    
    
    #conf_beta_global=np.linalg.pinv(conf_mat_global)@Xmat_deconf_global
    #conf_beta_global[np.abs(conf_beta_global)<1e-10]=0
    #Xmat_deconf_global=Xmat_deconf_global-conf_mat_global@conf_beta_global
    #Xmat_deconf_global=demean_mat(Xmat_deconf_global,0)
    Xmat=Xmat_pre.copy()   

    """
    phenotypes for supervision
    """    
    """
    cross-validated classification for targets
    """    
    for target_cnt,this_target in enumerate(dict_target_vectors.keys()):
        if target_cnt==required_target:
            
            print("********************************************")
            print("classification: "+this_target)
            print("********************************************")
            print("")
            target_type=copy.deepcopy(dict_target_types[this_target])
            cv_option=copy.deepcopy(cv_options[this_target])
            gaussianise_xmat=copy.deepcopy(gaussianise_dists[this_target][0])
            gaussianise_targ=copy.deepcopy(gaussianise_dists[this_target][1])
            gaussianise_cca=copy.deepcopy(gaussianise_dists[this_target][2])
            # if gaussianise_targ==False:#no target gaussianisation
            #     name_prefix=name_prefix+"_notargGaus"
            # if target_type in ['binary','multinomial','count']:
            #     target_deconf=False#if target is binary or multinomial don't deconfound it
            if target_type!='continuous':
                gaussianise_targ=False#if target is not continuous don't gaussianise it
            # target_deconf=False#if target is binary or multinomial don't deconfound it
            target_vector=copy.deepcopy(dict_target_vectors[this_target])
            nonan_positions=np.where(~np.isnan(target_vector))[0]
            target_vector_nonan=target_vector[nonan_positions].copy()#target vector may have a bunch of nan *subjects* (missing data), so we need to exclude them in xvalidation
            predicted_target=list()
            actual_target=list()
            actual_target_nodeconf=list()
            predicted_target_deconf_gx=list()
            predicted_target_deconf_lx=list()
            
            Xmat_nonan=Xmat[nonan_positions,:].copy()
            # if do_CCA:
            #     idpsup_mat_nonan=idp_supervision[nonan_positions,:].copy()#anything you do for Xmat_nonan, do for idpsup_mat_nonan too
            #     # n_cca=np.min([idp_supervision.shape[1],Xmat_nonan.shape[1]])
            #     ca = CCA(n_components=n_cca)
            
            if global_deconf:#double check this later
                print("prepare matrices for out-of-cross-validation deconfounding")
                """
                NOTE! Xmat is not demeaned, but Xmat_deconf_global and whatever goes into creating it is demeaned.
                """
                if gaussianise_xmat:                
                    conf_mat_global=conf_mat_raw.copy()#gaussianise_nonparametric(conf_mat_raw)
                    conf_mat_global=demean_mat(conf_mat_global,0)
                    
                    Xmat_nonan_normal_global,notofinterest=gaussianise_nonparametric(Xmat_nonan)
                    Xmat_nonan_deconf_global=demean_mat(Xmat_nonan_normal_global,0)

                    # if do_CCA:
                    #     idpsup_mat_nonan_normal_global,notofinterest=gaussianise_nonparametric(idpsup_mat_nonan)
                    #     idpsup_mat_nonan_deconf_global=demean_mat(idpsup_mat_nonan_normal_global,0)

                else:
                    conf_mat_global=demean_mat(conf_mat_raw,0)
                    Xmat_nonan_deconf_global=demean_mat(Xmat_nonan,0)
                    # if do_CCA:
                    #     idpsup_mat_nonan_deconf_global=demean_mat(idpsup_mat_nonan,0)
                
                if gaussianise_targ:
                    target_global_train_normal,target_trans_qt=gaussianise_nonparametric(target_vector_nonan[:,np.newaxis],nquantiles=int(len(target_vector_nonan)/(kvnum+1)))
                    
                Xmat_nonan_deconf_global,noint=deconfound_mat_fslnets(Xmat_nonan_deconf_global,conf_mat_global,train_beta=None)
                # if do_CCA:
                #     idpsup_mat_nonan_deconf_global,noint=deconfound_mat_fslnets(idpsup_mat_nonan_deconf_global,conf_mat_global,train_beta=None)
            
            for gb_kvcnt in global_kvnums:#range(global_kvnum):
                print("global kvnum: " + str(gb_kvcnt))
                randorder_subs=np.random.permutation(Xmat_nonan.shape[0])
                subnum_FoldPercent=int(len(randorder_subs)/kvnum)
                for kvcnt in range(kvnum):
                    print("kvnum: " +str(kvcnt))
                    test_subs=randorder_subs[kvcnt*subnum_FoldPercent:(kvcnt+1)*subnum_FoldPercent].copy()
                    train_subs = np.asarray([aa for aa in randorder_subs if aa not in test_subs])
                    
                    if xvalidated_featselect:
                        print("cross-validated feature selection based on correlation with target") 
                        corr2target=[np.corrcoef(Xmat_nonan[train_subs,ii],target_vector_nonan[train_subs])[0,1] for ii in range(Xmat_nonan.shape[1])]
                        if xvalidated_featselect_perc is not None:
                            argsort_corr2target=np.argsort(np.abs(corr2target))[::-1][:int(xvalidated_featselect_perc*len(corr2target))]
                        else:
                            argsort_corr2target=np.argsort(np.abs(corr2target))[::-1][:xvalidated_featselect_num]
                        Xmat_local=Xmat_nonan[:,argsort_corr2target].copy()  
                        print(Xmat_local.shape)
                    else:
                        Xmat_local=Xmat_nonan.copy()  
                    # if do_CCA:
                    #     if xvalidated_featselect:
                    #         print("cross-validated feature selection based on correlation with target") 
                    #         corr2target=[np.corrcoef(idpsup_mat_nonan[train_subs,ii],target_vector_nonan[train_subs])[0,1] for ii in range(idpsup_mat_nonan.shape[1])]
                    #         if xvalidated_featselect_perc is not None:
                    #             argsort_corr2target=np.argsort(np.abs(corr2target))[::-1][:int(xvalidated_featselect_perc*len(corr2target))]
                    #         else:
                    #             argsort_corr2target=np.argsort(np.abs(corr2target))[::-1][:xvalidated_featselect_num]
                    #         idpsup_mat_local=idpsup_mat_nonan[:,argsort_corr2target].copy()  
                    #         print(idpsup_mat_local.shape)
                    #     else:
                    #         idpsup_mat_local=idpsup_mat_nonan.copy()          
                    
                    if gaussianise_xmat:
                        print("gaussianise distributions: fold-wise")
                        #data
                        Xmat_local_train_normal,X_trans_qt=gaussianise_nonparametric(Xmat_local[train_subs,:],nquantiles=int(len(train_subs)/(kvnum+1)))
                        Xmat_local_test_normal,notofinterest=gaussianise_nonparametric(Xmat_local[test_subs,:],nquantiles=int(len(train_subs)/(kvnum+1)),X_trans_qt=X_trans_qt)
                    else:
                        Xmat_local_train_normal=Xmat_local[train_subs,:].copy()
                        Xmat_local_test_normal=Xmat_local[test_subs,:].copy()
                    # if do_CCA:
                    #     imp = IterativeImputer(max_iter=10, random_state=0)
                    #     imp.fit(idpsup_mat_local[train_subs,:])
                    #     # X_test = [[np.nan, 2], [6, np.nan], [np.nan, 6]]
                    #     idptrain_local=imp.transform(idpsup_mat_local[train_subs,:])
                    #     idptest_local=imp.transform(idpsup_mat_local[test_subs,:])
                    #     if gaussianise_cca:
                    #         print("gaussianise distributions: fold-wise")
                    #         #data
                    #         idpsup_local_train_normal,X_trans_qt=gaussianise_nonparametric(idptrain_local,nquantiles=int(len(train_subs)/(kvnum+1)))
                    #         idpsup_local_test_normal,notofinterest=gaussianise_nonparametric(idptest_local,nquantiles=int(len(train_subs)/(kvnum+1)),X_trans_qt=X_trans_qt)
                    #     else:
                    #         idpsup_local_train_normal=idptrain_local.copy()
                    #         idpsup_local_test_normal=idptest_local.copy()
                    if gaussianise_targ:
                        #target
                        target_local_train_normal,target_trans_qt=gaussianise_nonparametric(target_vector_nonan[train_subs][:,np.newaxis],nquantiles=int(len(train_subs)/(kvnum+1)))
                        target_local_test_normal,notofinterest=gaussianise_nonparametric(target_vector_nonan[test_subs][:,np.newaxis],nquantiles=int(len(train_subs)/(kvnum+1)),X_trans_qt=target_trans_qt)
                        target_local_train_normal=np.squeeze(target_local_train_normal)
                        target_local_test_normal=np.squeeze(target_local_test_normal)
    #                    #confounds
    #                    confmat_local_train_normal,conf_trans_qt=gaussianise_nonparametric(conf_mat_raw[train_subs,:],nquantiles=int(len(train_subs)/(kvnum+1)))
    #                    confmat_local_test_normal,notofinterest=gaussianise_nonparametric(conf_mat_raw[test_subs,:],nquantiles=int(len(train_subs)/(kvnum+1)),X_trans_qt=conf_trans_qt)
                    else:
                        target_local_train_normal=target_vector_nonan[train_subs].copy()
                        target_local_test_normal=target_vector_nonan[test_subs].copy()
                    
                    target_local_normal=target_local_test_normal.copy()
                        
                    confmat_local_train=conf_mat_raw[train_subs,:].copy()
                    confmat_local_test=conf_mat_raw[test_subs,:].copy()
                           
                    if xmat_deconf:
                        print("x-validated deconfounding of data matrix")
                        Xmat_deconf_local_train,train_X_beta_local=deconfound_mat_fslnets(Xmat_local_train_normal,confmat_local_train,train_beta=None)
                        Xmat_deconf_local_test,notint=deconfound_mat_fslnets(Xmat_local_test_normal,confmat_local_test,train_beta=train_X_beta_local)
                        # if do_CCA:
                        #     idpsup_deconf_local_train,train_sub_beta_local=deconfound_mat_fslnets(idpsup_local_train_normal,confmat_local_train,train_beta=None)
                        #     idpsup_deconf_local_test,notint=deconfound_mat_fslnets(idpsup_local_test_normal,confmat_local_test,train_beta=train_sub_beta_local)
                        
                    else:
                        Xmat_deconf_local_train=Xmat_local_train_normal.copy()
                        Xmat_deconf_local_test=Xmat_local_test_normal.copy()
                        # if do_CCA:
                        #     idpsup_deconf_local_train=idpsup_local_train_normal.copy()
                        #     idpsup_deconf_local_test=idpsup_local_test_normal.copy()
                        
                    if target_deconf:
                        print("x-validated deconfounding of target matrix")
                        target_deconf_local_train,train_targ_beta_local=deconfound_mat_fslnets(target_local_train_normal[:,np.newaxis],confmat_local_train,train_beta=None)
                        target_deconf_local_test,notint=deconfound_mat_fslnets(target_local_test_normal[:,np.newaxis],confmat_local_test,train_beta=train_targ_beta_local)
                        target_deconf_local_train=np.squeeze(target_deconf_local_train)
                        target_deconf_local_test=np.squeeze(target_deconf_local_test)
                    else:
                        target_deconf_local_train=target_local_train_normal.copy()
                        target_deconf_local_test=target_local_test_normal.copy()

                    if nodeconf_predict:
                        print("prediction- no deconfounding")
                        # if do_CCA: 
                        #     ca.fit(Xmat_local_train_normal, idpsup_local_train_normal)
                        #     Xmat_local_test_normal, idpsup_local_test_normal = ca.transform(Xmat_local_test_normal, idpsup_local_test_normal)
                        #     Xmat_local_train_normal, idpsup_local_train_normal = ca.transform(Xmat_local_train_normal, idpsup_local_train_normal)
                        print(Xmat_local_train_normal.shape)
                        print(Xmat_local_test_normal.shape)
                        print(target_local_train_normal.shape)
                        predtarget=cross_validated_regression(Xmat_local_train_normal,Xmat_local_test_normal,target_local_train_normal,target_type=target_type,options=cv_option)#,kvnormalize=kvnormalize)
                        predicted_target.append(predtarget.copy())    
                        actual_target_nodeconf.append(target_local_normal)                   
                    
                    if global_deconf:
                        print("prediction- global deconfounding X")#CCA doest work with this
                        predtarget_deconf_gx=cross_validated_regression(Xmat_nonan_deconf_global[train_subs,:],Xmat_nonan_deconf_global[test_subs,:],target_global_train_normal[train_subs],target_type=target_type,options=cv_option)#,kvnormalize=kvnormalize)
                        predicted_target_deconf_gx.append(predtarget_deconf_gx.copy())
                    
                    print("prediction- local deconfounding X")
                    # if do_CCA: 
                    #     ca.fit(Xmat_deconf_local_train, idpsup_deconf_local_train)
                    #     Xmat_deconf_local_test, idpsup_deconf_local_test = ca.transform(Xmat_deconf_local_test, idpsup_deconf_local_test)
                    #     Xmat_deconf_local_train, idpsup_deconf_local_train = ca.transform(Xmat_deconf_local_train, idpsup_deconf_local_train)
                    predtarget_deconf_lx=cross_validated_regression(Xmat_deconf_local_train,Xmat_deconf_local_test,target_deconf_local_train,target_type=target_type,options=cv_option)#,kvnormalize=kvnormalize)
                    predicted_target_deconf_lx.append(predtarget_deconf_lx.copy())    
                    
                    if target_type=='binary':
                        thisvar=target_deconf_local_test.copy()                        
                        thisvar[thisvar<np.mean(thisvar)]=0
                        thisvar[thisvar>0]=1
                        actual_target.append(thisvar.copy())
                    else:
                        actual_target.append(target_deconf_local_test.copy())
                                        
                #save
                if nodeconf_predict:
                    outname=thisclass_path_stoch+name_prefix+"_predicted_"+this_target+"_"+target_type+"_"+str(gb_kvcnt)+"_"+str(global_kvnum)+".npy"#gb_kvcnt
                    np.save(outname,predicted_target)

                    outname=thisclass_path_stoch+name_prefix+"_actual_nodeconf_"+this_target+"_"+target_type+"_"+str(gb_kvcnt)+"_"+str(global_kvnum)+".npy"
                    np.save(outname,actual_target_nodeconf)
                
                outname=thisclass_path_stoch+name_prefix+"_actual_"+this_target+"_"+target_type+"_"+str(gb_kvcnt)+"_"+str(global_kvnum)+".npy"
                np.save(outname,actual_target)
                print(outname)
                
                if global_deconf:
                    outname=thisclass_path_stoch+name_prefix+"_predicted_"+this_target+"_"+target_type+"_deconf_gx_"+str(gb_kvcnt)+"_"+str(global_kvnum)+".npy"
                    np.save(outname,predicted_target_deconf_gx)
                
                outname=thisclass_path_stoch+name_prefix+"_predicted_"+this_target+"_"+target_type+"_deconf_lx_"+str(gb_kvcnt)+"_"+str(global_kvnum)+".npy"
                np.save(outname,predicted_target_deconf_lx)