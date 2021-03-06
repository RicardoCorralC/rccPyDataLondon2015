# Ricardo Corral Corral, June 2015
# Easy model evaluation over CATH and SCOP datasets

import sys
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.cross_validation import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn import grid_search
from collections import Counter
import matplotlib.pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from sklearn.grid_search import RandomizedSearchCV

def dbParser(db,level=1,dbtype='CATH',minsamples=500,verbose=False):
    sepchar = ''
    if dbtype == 'CATH':
        sepchar = '_'
    elif dbtype == 'SCOP':
        sepchar = '.'

    _c = Counter()
    domnames, vects, classifs = [], [], []
    domnames2, vects2, classifs2 = [], [], []

    fin = open(db)
    for l in fin:
        l = l.strip().split(',')
        domnames.append(l[0])
        vects.append( map(float,l[1:-1]) )
        _class = '_'.join(l[-1].split( sepchar )[:level])
        classifs.append(  _class  )
        _c[_class] += 1

    skippedclasses = set()
    for i in xrange(len(domnames)):
        if _c[classifs[i]] < minsamples:
            skippedclasses.add( classifs[i] )
            continue
        domnames2.append( domnames[i] )
        vects2.append( vects[i] )
        classifs2.append( classifs[i] )

    if verbose:
        print '\nClassess skipped having less than ', minsamples, 'samples:\n', '\n'.join(skippedclasses) , '\n'


    dataDict = dict()
    dataDict['domain_names'] = np.asarray(domnames2)
    dataDict['vectors'] = np.asarray(vects2)
    dataDict['target_names'] = np.asarray(classifs2)
    return dataDict



def makeElenaTuning(dataset,dbtype='CATH',level=1,k_iters=3,clf = ExtraTreesClassifier(n_estimators=5,class_weight='auto')):
    #Adaptation from https://github.com/cambridgecoding/pydata-tutorial
    dataDict = dbParser(dataset,level=level,dbtype=dbtype)
    print dataDict

    labels = dataDict['target_names']

    parameters = {"min_samples_split":[1,2,3,4,5], "max_features":[None,1,2,3,4,5]}


    grid = grid_search.GridSearchCV(clf, parameters, cv= 3)
    grid.fit(dataDict['vectors'], dataDict['target_names'])

    print "The best parameters are: min_samples_split=", grid.best_params_['min_samples_split'],"and max_depth=",grid.best_params_['max_features']

    score_dict = grid.grid_scores_
    scores = [x[1] for x in score_dict]
    scores = np.array(scores).reshape(len(parameters["min_samples_split"]), len(parameters["max_features"]))
    scores = np.transpose(scores)

    # Make a heatmap with the performance
    #plt.figure(figsize=(12, 6))
    plt.imshow(scores, interpolation='nearest', origin='higher', cmap=plt.cm.get_cmap('jet_r'))
    plt.xticks(np.arange(len(parameters["min_samples_split"])), parameters["min_samples_split"])
    plt.yticks(np.arange(len(parameters["max_features"])), parameters["max_features"])
    plt.xlabel('min_samples_split')
    plt.ylabel('max_features')

    cbar = plt.colorbar()
    cbar.set_label('Classification Accuracy', rotation=270, labelpad=20)

    plt.show()


def report(grid_scores, n_top=3):
    from operator import itemgetter
    top_scores = sorted(grid_scores, key=itemgetter(1), reverse=True)[:n_top]
    for i, score in enumerate(top_scores):
        print("Model with rank: {0}".format(i + 1))
        print("Mean validation score: {0:.3f} (std: {1:.3f})".format(
              score.mean_validation_score,
              np.std(score.cv_validation_scores)))
        print("Parameters: {0}".format(score.parameters))
        print("")


def makeRandomCV(dataset,dbtype='CATH',
                level=1,
                k_iters=10,
                minsamples=500,
                clf = ExtraTreesClassifier(n_estimators=5,class_weight='auto')):

    from scipy.stats import randint as sp_randint
    
    dataDict = dbParser(dataset,level=level,dbtype=dbtype,minsamples=minsamples)
    print dataDict

    labels = dataDict['target_names']

    param_dist = {"max_depth": [3, None],
              "max_features": sp_randint(1, 11),
              "min_samples_split": sp_randint(1, 11),
              "min_samples_leaf": sp_randint(1, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

    n_iter_search = k_iters
    random_search = RandomizedSearchCV(clf, param_distributions=param_dist,
                                   n_iter=n_iter_search)


    random_search.fit(dataDict['vectors'], dataDict['target_names'])
    report(random_search.grid_scores_)




def makeAllEvals(dataset,dbtype='CATH',level=1,k_iters=10):

    dataDict = dbParser(dataset,level=level,dbtype=dbtype)
    print dataDict

    labels = dataDict['target_names']
    skf = StratifiedKFold(labels, k_iters)


    clf = ExtraTreesClassifier(n_estimators=250 ,class_weight='auto',max_depth=None,min_samples_split=4)

    accsList = []

    for train, test in skf:
        print '\n--------------------------------------------------\n'
        _train = [dataDict['vectors'][i] for i in train]
        _test = [dataDict['vectors'][i] for i in test]
        _targets = [dataDict['target_names'][i] for i in train]

        clf.fit(_train,_targets)

        y_true = [labels[i] for i in test]
        y_pred = clf.predict(_test)


        localAccuracy = accuracy_score(y_true, y_pred)
        accsList.append(localAccuracy)
        print '*** localAccuracy:', localAccuracy , '\n'

        print(classification_report(y_true, y_pred))

        cm = confusion_matrix(y_true, y_pred,labels=clf.classes_)
        print cm

    _ACC = np.mean(accsList)
    _ACCstd = np.std(accsList)

    print '\n[ FINAL SUMMARY ]'
    print ' *** ACCURACY: ', _ACC
    print ' *** ACCURACY deviation: ', _ACCstd



def main():
    fn = sys.argv[1] #textfile with each line being: domainname, 26dvector, classification
    level = int(sys.argv[2]) # An integer number
    dbtype = sys.argv[3] # CATH or SCOP

    makeElenaTuning(dataset=fn,level=level,dbtype=dbtype)
    exit()
    makeAllEvals(dataset=fn,level=level,dbtype=dbtype)

if __name__ == '__main__':
    main()
