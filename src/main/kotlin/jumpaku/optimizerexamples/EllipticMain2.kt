package jumpaku.optimizerexamples

import jumpaku.commons.control.result
import jumpaku.curves.core.curve.Curve
import jumpaku.curves.core.curve.Interval
import jumpaku.curves.core.curve.arclength.ReparametrizedCurve
import jumpaku.curves.core.curve.bspline.BSpline
import jumpaku.curves.core.geom.Point
import jumpaku.curves.core.geom.lerp
import jumpaku.curves.core.transform.Calibrate
import kotlin.math.acos
import kotlin.math.cos
import kotlin.math.sin


fun evaluateFarsAndWeight2(xs: DoubleArray, s: BSpline, farsCount: Int): EvalResult {
    require(xs.size == 4)
    val weight = (sigmoid(xs[0]) * 2 - 1).coerceIn(-0.999, 0.999)
    val theta = 2 * acos(weight)

    val circularArc = object : Curve {
        override val domain: Interval = Interval(0.0, 2 * theta)
        override fun evaluate(t: Double): Point = Point.xy(cos(t), sin(t))
    }
    val ps = circularArc.evaluateAll(5)

    val (t0, t4) = s.domain
    val t2 = t0.lerp(sigmoid(xs[1]), t4).coerceIn(t0, t4)
    val t1 = t0.lerp(sigmoid(xs[2]), t2).coerceIn(t0, t2)
    val t3 = t2.lerp(sigmoid(xs[3]), t4).coerceIn(t2, t4)
    val ts = listOf(t0, t1, t2, t3, t4)
    println(ts)
    val qs = ts.map(s)

    val calibrate = result {
        Calibrate(ps[1] to qs[1], ps[2] to qs[2], ps[3] to qs[3])
    }.tryRecover {
        Calibrate(ps[1] to qs[1], ps[3] to qs[3])
    }.orRecover {
        Calibrate(ps[2] to qs[2])
    }

    val transformed = object : Curve {
        override val domain: Interval = circularArc.domain
        override fun evaluate(t: Double): Point = calibrate(circularArc(t))
    }
    val samples = 100
    val reparametrizedCA = ReparametrizedCurve.of(transformed, transformed.domain.sample(samples))
    val reparametrisedFSC = ReparametrizedCurve.of(s, s.domain.sample(samples))
    val ls = transformed.domain.sample(farsCount).map { reparametrizedCA.reparametrizer.toArcLengthRatio(it) }
    val rCA = ls.map(reparametrizedCA)
    val rFSC = ls.map(reparametrisedFSC)
    val error = rCA.zip(rFSC).sumByDouble { (c, f) -> c.distSquare(f) }
    val params = ls.map { reparametrisedFSC.reparametrizer.toOriginal(it) }
    return EvalResult(weight, params, error)
}
