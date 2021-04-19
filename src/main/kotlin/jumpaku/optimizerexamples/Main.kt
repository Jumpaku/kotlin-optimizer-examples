package jumpaku.optimizerexamples

import org.apache.commons.math3.linear.MatrixUtils
import org.apache.commons.math3.linear.RealVector
import org.apache.commons.math3.optim.InitialGuess
import org.apache.commons.math3.optim.MaxEval
import org.apache.commons.math3.optim.SimpleBounds
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer
import org.apache.commons.math3.random.MersenneTwister
import kotlin.math.*


fun cmaes(seed: Int, iterations: Int, initial: RealVector, objective: (RealVector) -> Double): Pair<RealVector, Double> {
    val dimension = initial.dimension
    val populationSize = 5* (4 + 3 * ceil(log(E, dimension.toDouble()))).toInt()
    val optimizer = CMAESOptimizer(
            iterations,
            0.0,
            true,
            0,
            0,
            MersenneTwister(seed),
            true,
            { iteration, _, _ -> iteration >= iterations }
    )
    val result = optimizer.optimize(
            CMAESOptimizer.Sigma(DoubleArray(dimension) { 1.0 }),
            CMAESOptimizer.PopulationSize(populationSize),
            ObjectiveFunction { x -> objective(MatrixUtils.createRealVector(x)) },
            GoalType.MINIMIZE,
            InitialGuess(initial.toArray()),
            SimpleBounds(
                    DoubleArray(dimension) { Double.NEGATIVE_INFINITY },
                    DoubleArray(dimension) { Double.POSITIVE_INFINITY }),
            MaxEval(populationSize * iterations)
    )
    optimizer.statisticsFitnessHistory.forEachIndexed { i, f ->
        println("$i: $f")
    }
    return result.run { MatrixUtils.createRealVector(point) to value }
}

fun rastrigin(x: RealVector): Double {
    val n = x.dimension
    val a = 10
    val fx = a * n + x.toArray().sumByDouble { xi -> xi * xi - a * cos(2 * PI * xi) }
    return fx
}

fun ackley(x: RealVector): Double {
    val a = -20 * exp(-0.2 * sqrt(x.dotProduct(x) / 2))
    val b = -exp(0.5 * cos(2 * PI * x.getEntry(0)) + 0.5 * cos(2 * PI * x.getEntry(1)))
    val fx = a + b + E + 20
    return fx
}

fun rosenbrock(x: RealVector): Double {
    val fx = x.toArray().toList().zipWithNext { a, b -> 100 * (b - a * a) * (b - a * a) + (1 - a) * (1 - a) }.sum()
    return fx
}

fun main() {
    val initial = MatrixUtils.createRealVector(doubleArrayOf(-1.0, -1.0))
    val answer = cmaes(
            1089,
            50,
            initial,
            ::rastrigin
    )
    println(answer)
}