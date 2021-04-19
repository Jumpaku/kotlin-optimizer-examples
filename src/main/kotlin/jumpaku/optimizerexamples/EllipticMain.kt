package jumpaku.optimizerexamples

import jumpaku.commons.control.None
import jumpaku.commons.control.Option
import jumpaku.commons.control.Some
import jumpaku.commons.control.result
import jumpaku.curves.core.curve.Interval
import jumpaku.curves.core.curve.ParamPoint
import jumpaku.curves.core.curve.arclength.ReparametrizedCurve
import jumpaku.curves.core.curve.bspline.BSpline
import jumpaku.curves.core.geom.Point
import jumpaku.curves.core.geom.lerp
import jumpaku.curves.core.transform.Calibrate
import jumpaku.curves.core.transform.Translate
import jumpaku.curves.fsc.DrawingStroke
import jumpaku.curves.fsc.generate.Fuzzifier
import jumpaku.curves.fsc.generate.Generator
import jumpaku.curves.fsc.identify.primitive.reference.EllipticGenerator
import jumpaku.curves.graphics.drawCubicBSpline
import jumpaku.curves.graphics.drawPoints
import jumpaku.curves.graphics.swing.DrawingPanel
import org.apache.commons.math3.optim.*
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer
import org.apache.commons.math3.random.MersenneTwister
import org.apache.commons.math3.util.ArithmeticUtils
import java.awt.Color
import java.awt.Dimension
import java.awt.Graphics
import java.awt.Graphics2D
import javax.swing.JFrame
import javax.swing.JPanel
import javax.swing.SwingUtilities
import kotlin.math.*

fun main() {
    SwingUtilities.invokeLater {
        val demo = DemoPanel()
        val drawing = DrawingPanel().apply {
            addCurveListener {
                val w = -0.707//0.999//-0.5//0.0//0.707//sqrt(0.5+0.5*0.707)
                val theta = 2 * 2 * acos(w)
                val span = 10.0
                val samples = 200
                val r = 100.0
                val stroke = Interval(PI / 2 - theta / 2, PI / 2 + theta / 2).sample(samples).mapIndexed { index, a ->
                    val param = span * index / samples
                    val point = Point.xy(
                            //r * cos(a), -r * sin(a)
                            param * 20, 0.0
                    ).transform(Translate(r * 2, r * 2))
                    ParamPoint(point, param)
                }
                //demo.update(DrawingStroke(stroke))
                demo.update(it.drawingStroke)
            }
            add(demo)
        }
        JFrame("GestureDemo").apply {
            defaultCloseOperation = JFrame.EXIT_ON_CLOSE
            contentPane.add(drawing)
            pack()
            isVisible = true
        }
    }
}

object Settings {

    val width = 640

    val height = 480

    val generator: Generator = Generator(
            degree = 3,
            knotSpan = 0.1,
            fillSpan = 0.0375,
            extendInnerSpan = 0.075,
            extendOuterSpan = 0.075,
            extendDegree = 2,
            fuzzifier = Fuzzifier.Linear(
                    velocityCoefficient = 0.007,
                    accelerationCoefficient = 0.008
            )
    )
}

class DemoPanel : JPanel() {

    init {
        preferredSize = Dimension(
                Settings.width,
                Settings.height
        )
    }

    var fscOpt: Option<BSpline> = None
    var updated = false

    fun update(drawingStroke: DrawingStroke) {
        val fsc = Settings.generator.generate(drawingStroke)

        fscOpt = Some(fsc)
        updated = true

        repaint()
    }

    override fun paint(g: Graphics) {
        if (!updated) return
        fscOpt.forEach { s ->
            val (w, ts) = findFarsAndWeight(0, s)
            val ps = ts.map(s)
            println("params: $ts")
            println("weight: $w")
            val g2d = g as Graphics2D
            g2d.drawCubicBSpline(s) { it.color = Color.BLACK }
            g2d.drawPoints(ts.map(s).map { it.copy(r = 5.0) }) { it.color = Color.BLACK }
            g2d.drawPoints(listOf(
                    bezier(-0.5 / w, ps[1], ps[2], ps[3], w).copy(r = 6.0),
                    bezier((w + 0.5) / w, ps[1], ps[2], ps[3], w).copy(r = 6.0)
            ))
            updated = false
        }
    }
}

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

fun evaluateFarsAndWeight(xs: DoubleArray, s: BSpline, farsCount: Int): EvalResult {
    require(xs.size == 4)
    val weight = (sigmoid(xs[0]) * 2 - 1).coerceIn(-0.999, 0.999)
    val (t0, t4) = s.domain
    val t2 = t0.lerp(sigmoid(xs[1]), t4).coerceIn(t0, t4)
    val t1 = t0.lerp(sigmoid(xs[2]), t2).coerceIn(t0, t2)
    val t3 = t2.lerp(sigmoid(xs[3]), t4).coerceIn(t2, t4)
    val ts = listOf(t0, t1, t2, t3, t4)
    val ps = ts.map(s)
    val theta = 2 * acos(weight)
    val qs = (0 until farsCount).map { i ->
        val a = i * theta * 0.5 + PI / 2
        Point.xy(cos(a), sin(a))
    }
    val calibrate = result {
        Calibrate(ps[1] to qs[1], ps[2] to qs[2], ps[3] to qs[3])
    }.tryRecover {
        Calibrate(ps[1] to qs[1], ps[3] to qs[3])
    }.orRecover {
        Calibrate(ps[2] to qs[2])
    }
    val reparametrised = ReparametrizedCurve.of(s.transform(calibrate), s.domain.sample(100))
    val rs = reparametrised.sample(farsCount)
    val params = rs.map { reparametrised.reparametrizer.toOriginal(it.param) }
    val error = qs.zip(rs).sumByDouble { (q, r) -> q.distSquare(r.point) }
    return EvalResult(weight, params, error)
}

fun findFarsAndWeight(depth: Int, s: BSpline): Pair<Double, List<Double>> {
    val nSamples = 100
    val (t0, t4) = s.domain
    val t2 = EllipticGenerator.computeEllipticFar(s, t0, t4, nSamples)
    val t1 = EllipticGenerator.computeEllipticFar(s, t0, t2, nSamples)
    val t3 = EllipticGenerator.computeEllipticFar(s, t2, t4, nSamples)
    val w = EllipticGenerator.computeEllipticWeight(s, t1, t3, t2, s.domain, nSamples)
    val initial = doubleArrayOf(
            log((1 + w) / (1 - w), E),
            log(((t2 - t0) / (t4 - t0)).let { y -> y / (1 - y) }, E),
            log(((t1 - t0) / (t2 - t0)).let { y -> y / (1 - y) }, E),
            log(((t3 - t2) / (t4 - t2)).let { y -> y / (1 - y) }, E)
    )
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
    val result = optimizer.optimize(
            CMAESOptimizer.Sigma(doubleArrayOf(0.01, 0.01, 0.01, 0.01)),
            CMAESOptimizer.PopulationSize(populationSize),
            ObjectiveFunction { x -> evaluateFarsAndWeight(x, s, farsCount).evaluation },
            GoalType.MINIMIZE,
            InitialGuess(initial),
            SimpleBounds.unbounded(dimension),
            MaxEval(populationSize * iterations + 1)
    ).run { evaluateFarsAndWeight(point, s, farsCount) }
    println("$w, ${result.weight}")
    (optimizer.statisticsFitnessHistory.forEach(::println))
    return result.run { this.weight to params }
}