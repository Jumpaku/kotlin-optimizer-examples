package jumpaku.optimizerexamples

import jumpaku.curves.core.curve.bspline.BSpline
import jumpaku.curves.core.geom.Point
import jumpaku.curves.core.geom.lerp
import jumpaku.curves.fsc.identify.primitive.reference.EllipticGenerator
import org.apache.commons.math3.optim.*
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer
import org.apache.commons.math3.random.MersenneTwister
import org.apache.commons.math3.util.ArithmeticUtils
import kotlin.math.E
import kotlin.math.exp
import kotlin.math.log
import kotlin.math.sqrt


data class EvalResult(val weight: Double, val params: List<Double>, val evaluation: Double)

fun sigmoid(x: Double): Double = (1 / (1 + exp(-x)))
fun bezier(t: Double, p0: Point, p1: Point, p2: Point, w: Double): Point {
    val c = 1 / (1 - w)
    if (c.isFinite() && t !in 0.0..1.0) {
        val q1 = p1.lerp(c to p0, c to p2)
        return bezier(1 / (2 - 1 / t), p0, q1, p2, -w)
    }
    val wt = ((1 - t) * (1 - t)) + (2 * t * (1 - t) * w) + (t * t)
    val c0 = ((1 - t) * (1 - t) - t * (1 - t)) / wt
    val c2 = (t * t - t * (1 - t)) / wt
    return p1.lerp(c0 to p0, c2 to p2)
}

fun findFarsAndWeight(depth: Int, s: BSpline): Pair<Double, List<Double>> {
    val nSamples = 100
    val (t0, t4) = s.domain
    val t2 = EllipticGenerator.computeEllipticFar(s, t0, t4, nSamples)
    val t1 = EllipticGenerator.computeEllipticFar(s, t0, t2, nSamples)
    val t3 = EllipticGenerator.computeEllipticFar(s, t2, t4, nSamples)
    val w = EllipticGenerator.computeEllipticWeight(s, t1, t3, t2, s.domain, nSamples)
    val initial = doubleArrayOf(0.0, 0.0, 0.0, 0.5 * sqrt(2.0))
    /*doubleArrayOf (
            log((1 + w) / (1 - w), E),
    log(((t2 - t0) / (t4 - t0)).let { y -> y / (1 - y) }, E),
    log(((t1 - t0) / (t2 - t0)).let { y -> y / (1 - y) }, E),
    log(((t3 - t2) / (t4 - t2)).let { y -> y / (1 - y) }, E)
    )*/
    val iterations = 100
    val checker = ConvergenceChecker<PointValuePair> { iteration, _, _ -> iteration >= iterations }
    val optimizer = CMAESOptimizer(
        iterations,
        0.0,
        false,
        0,
        1,
        MersenneTwister(1089),
        true,
        checker
    )
    val weight = (sigmoid(initial[0]) * 2 - 1).coerceIn(-0.999, 0.999)
    println("$weight, $w")
    val (s0, s4) = s.domain
    val s2 = t2
    val s1 = s0.lerp(sigmoid(initial[1]), s2).coerceIn(s0, s2)
    val s3 = s2.lerp(sigmoid(initial[2]), s4).coerceIn(s2, s4)
    println(listOf(t0, t1, t2, t3, t4))
    println(listOf(s0, s1, s2, s3, s4))
    val farsCount = 1 + ArithmeticUtils.pow(2, depth + 2)
    val dimension = 4
    val populationSize = 100
    val r = optimizer.optimize(
        CMAESOptimizer.Sigma(doubleArrayOf(1.0, 1.0, 1.0, 1.0)),
        CMAESOptimizer.PopulationSize(populationSize),
        ObjectiveFunction { x -> evaluateFarsAndWeight2(x, s, farsCount).evaluation },
        GoalType.MINIMIZE,
        InitialGuess(initial),
        SimpleBounds.unbounded(dimension),
        MaxEval(populationSize * iterations + 1)
    )
    val result = r.run { evaluateFarsAndWeight2(point, s, farsCount) }
    println("$w, ${result.weight}, ${result.evaluation}")

    (optimizer.statisticsFitnessHistory.forEach(::println))
    return result.run { this.weight to params }
}