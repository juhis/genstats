'use strict'

var descriptives = require('./descriptives')
var ttest = require('./ttest')
var wilcoxon = require('./wilcoxon')

module.exports = {

    mean: descriptives.mean,
    variance: descriptives.variance,
    stdev: descriptives.stdev,

    covariance: descriptives.covariance,
    correlation: descriptives.correlation,
    correlationStdNorm: descriptives.correlationStdNorm,

    standardNormalize: descriptives.standardNormalize,
    
    welch: ttest.welch,
    student: ttest.student,
    wilcoxon: wilcoxon.wilcoxon,
    wilcoxonRanks: wilcoxon.wilcoxonRanks
    
}
