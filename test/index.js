var mocha = require('mocha')
var assert = require('chai').assert
var genstats = require('../index')

var arr = [-1, -0.5, 0]

describe('genstats#mean', function() {

    it('returns -0.5 for [-1, -0.5, 0]', function() {
        assert.equal(-0.5, genstats.mean(arr))
    })
})

describe('genstats#variance', function() {

    it('returns 0.25 for [-1, -0.5, 0]', function() {
        assert.equal(0.25, genstats.variance(arr))
    })
})

describe('genstats#stdev', function() {

    it('returns 0.5 for [-1, -0.5, 0]', function() {
        assert.equal(0.5, genstats.stdev(arr))
    })
})

describe('genstats#covariance', function() {

    it('returns 0.25 for [-1, -0.5, 0], [0, 0.5, 1]', function() {
        assert.equal(0.25, genstats.covariance(arr, [0, 0.5, 1]))
    })
})

describe('genstats#correlation', function() {

    it('returns -1 for [1, 0.5, 0], [-1, -0.5, 0]', function() {
        assert.equal(-1, genstats.correlation([1, 0.5, 0], arr))
    })
    it('returns about 0.69 for [-1, -0.5, 1], [-1, 0, 0]', function() {
        assert.isAbove(0.7, genstats.correlation([-1, -0.5, 1], [-1, 0, 0]))
        assert.isBelow(0.69, genstats.correlation([-1, -0.5, 1], [-1, 0, 0]))
    })
})

describe('genstats#correlationStdNorm', function() {

    it('returns -0.625 for [1, 0.5, 0], [0, -0.5, 1]', function() {
        assert.equal(-0.625, genstats.correlationStdNorm([1, 0.5, 0], arr))
    })
})

describe('genstats#standardNormalize', function() {

    it('turns [-1, -0.5, 0] to [-1, 0, 1]', function() {
        var a = arr.slice(0)
        genstats.standardNormalize(a)
        assert.deepEqual(a, [-1, 0, 1])
    })
})

describe('genstats#student', function() {

    it('returns {t: 0, p: 0.5, df: 4} for [-1, -0.5, 0], [-1, -0.5, 0]', function() {
        assert.deepEqual({t: 0, p: 0.5, df: 4}, genstats.student(arr, arr))
    })

    it('returns {-0.78 < t < -0.77, 0.240 < p < 0.241, df: 4} for [-1, -0.5, 0], [-0.5, -0.25, 0]', function() {
        var result = genstats.student(arr, [-0.5, -0.25, 0])
        assert.isAbove(-0.77, result.t)
        assert.isBelow(-0.78, result.t)
        assert.isAbove(0.241, result.p)
        assert.isBelow(0.240, result.p)
        assert.equal(4, result.df)
    })

    it('returns a 1e-100 < p < 1e-30 for wildly different random arrays of length 100', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 100; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.isAbove(1e-30, genstats.student(a1, a2).p)
        assert.isBelow(1e-100, genstats.student(a1, a2).p)
    })

    it('returns a p of Number.MIN_VALUE for wildly different random arrays of length 1000', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 1000; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.equal(Number.MIN_VALUE, genstats.student(a1, a2).p)
    })
})

describe('genstats#welch', function() {

    it('returns {t: 0, p: 0.5, df: 4} for [-1, -0.5, 0], [-1, -0.5, 0]', function() {
        assert.deepEqual({t: 0, p: 0.5, df: 4}, genstats.welch(arr, arr))
    })
    it('returns {-0.78 < t < -0.77, 0.248 < p < 0.249, 2.9 < df < 3} for [-1, -0.5, 0], [-0.5, -0.25, 0]', function() {
        var result = genstats.welch(arr, [-0.5, -0.25, 0])
        assert.isAbove(-0.77, result.t)
        assert.isBelow(-0.78, result.t)
        assert.isAbove(0.249, result.p)
        assert.isBelow(0.248, result.p)
        assert.isAbove(3, result.df)
        assert.isBelow(2.9, result.df)
    })

    it('returns a 1e-100 < p < 1e-30 for wildly different random arrays of length 100', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 100; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.isAbove(1e-30, genstats.welch(a1, a2).p)
        assert.isBelow(1e-100, genstats.welch(a1, a2).p)
    })

    it('returns a p of Number.MIN_VALUE for wildly different random arrays of length 1000', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 1000; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.equal(Number.MIN_VALUE, genstats.welch(a1, a2).p)
    })
})

describe('genstats#wilcoxon', function() {

    it('returns {0.82 < p < 0.83, auc: 1/3} for [-1, -0.5, 0], [-1, -0.5, 0]', function() {
        var result = genstats.wilcoxon(arr, arr)
        assert.isAbove(0.83, result.p)
        assert.isBelow(0.82, result.p)
        assert.equal(1/3, result.auc)
    })

    it('returns {p: 0.66 < p < 0.67, auc: 2/9} for [-1, -0.5, 0], [-0.5, -0.25, 0]', function() {
        var result = genstats.wilcoxon(arr, [-0.5, -0.25, 0])
        assert.isAbove(0.67, result.p)
        assert.isBelow(0.66, result.p)
        assert.equal(2/9, result.auc)
    })

    it('returns a 1e-100 < p < 1e-30 for wildly different random arrays of length 100', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 100; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.isAbove(1e-30, genstats.wilcoxon(a1, a2).p)
        assert.isBelow(1e-100, genstats.wilcoxon(a1, a2).p)
    })

    it('returns a p of Number.MIN_VALUE for wildly different random arrays of length 1000', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 1000; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.equal(Number.MIN_VALUE, genstats.wilcoxon(a1, a2).p)
    })
})

describe('genstats#wilcoxon', function() {

    it('returns {0.82 < p < 0.83, auc: 1/3} for [-1, -0.5, 0], [-1, -0.5, 0]', function() {
        var result = genstats.wilcoxon(arr, arr)
        assert.isAbove(0.83, result.p)
        assert.isBelow(0.82, result.p)
        assert.equal(1/3, result.auc)
    })
    it('returns {p: 0.66 < p < 0.67, auc: 2/9} for [-1, -0.5, 0], [-0.5, -0.25, 0]', function() {
        var result = genstats.wilcoxon(arr, [-0.5, -0.25, 0])
        assert.isAbove(0.67, result.p)
        assert.isBelow(0.66, result.p)
        assert.equal(2/9, result.auc)
    })

    it('returns a 1e-100 < p < 1e-30 for wildly different random arrays of length 100', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 100; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.isAbove(1e-30, genstats.wilcoxon(a1, a2).p)
        assert.isBelow(1e-100, genstats.wilcoxon(a1, a2).p)
    })

    it('returns a p of Number.MIN_VALUE for wildly different random arrays of length 1000', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 1000; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.equal(Number.MIN_VALUE, genstats.wilcoxon(a1, a2).p)
    })
})

describe('genstats#wilcoxonRanks', function() {

    it('returns {0.82 < p < 0.83, auc: 1/3} for [0, 1, 2], [0, 1, 2]', function() {
        var ranks1 = [0, 1, 2]
        var ranks2 = [0, 1, 2]
        var result = genstats.wilcoxon(ranks1, ranks2)
        assert.isAbove(0.83, result.p)
        assert.isBelow(0.82, result.p)
        assert.equal(1/3, result.auc)
    })

    it('returns {p: 0.66 < p < 0.67, auc: 2/9} for [0, 1, 3], [1, 2, 3]', function() {
        var ranks1 = [0, 1, 3]
        var ranks2 = [1, 2, 3]
        var result = genstats.wilcoxon(ranks1, ranks2)
        assert.isAbove(0.67, result.p)
        assert.isBelow(0.66, result.p)
        assert.equal(2/9, result.auc)
    })

    it('returns a 1e-100 < p < 1e-30 for wildly different random arrays of length 100', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 100; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.isAbove(1e-30, genstats.wilcoxon(a1, a2).p)
        assert.isBelow(1e-100, genstats.wilcoxon(a1, a2).p)
    })

    it('returns a p of Number.MIN_VALUE for wildly different random arrays of length 1000', function() {
        var a1 = []
        var a2 = []
        for (var i = 0; i < 1000; i++) {
            a1.push(-10000 * Math.random())
            a2.push(10000 * Math.random())
        }
        assert.equal(Number.MIN_VALUE, genstats.wilcoxon(a1, a2).p)
    })
})
