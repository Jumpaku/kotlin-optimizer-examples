package jumpaku.optimizerexamples

import jumpaku.commons.control.None
import jumpaku.commons.control.Option
import jumpaku.commons.control.Some
import jumpaku.curves.core.curve.Interval
import jumpaku.curves.core.curve.ParamPoint
import jumpaku.curves.core.curve.bspline.BSpline
import jumpaku.curves.core.geom.Point
import jumpaku.curves.core.transform.Translate
import jumpaku.curves.fsc.DrawingStroke
import jumpaku.curves.fsc.generate.Fuzzifier
import jumpaku.curves.fsc.generate.Generator
import jumpaku.curves.graphics.DrawStyle
import jumpaku.curves.graphics.drawCubicBSpline
import jumpaku.curves.graphics.drawPoints
import jumpaku.curves.graphics.swing.DrawingPanel
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
                val w = 0.1//sqrt(0.5+0.5*0.707)//0.707//0.0//0.707//
                val theta = 2 * 2 * acos(w)
                val span = 10.0
                val samples = 200
                val r = 100.0
                val stroke = Interval(PI / 2 - theta / 2, PI / 2 + theta / 2).sample(samples).mapIndexed { index, a ->
                    val param = span * index / samples
                    val point = Point.xy(
                        r * cos(a), -r * sin(a)
                        //param * 20, 0.0
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
            g2d.drawPoints(ts.map(s).mapIndexed { i, p -> p.copy(r = 2.0 * (i + 1)) }) { it.color = Color.BLACK }
            val evals = Interval(min(0.0, -0.5 / w), max(1.0, (w + 0.5) / w))
                .sample(0.05)
                .map { bezier(it, ps[1], ps[2], ps[3], w) }
            g2d.drawPoints(evals.map { it.copy(r = 2.0) }, DrawStyle(Color.RED))
            updated = false
        }
    }
}