using System;

for (int i = 1; i < 5; i++)
{
    for (int j = i; j < 5; j++)
    {
        for (double u = 3; u <= 6; u += 3)
        {
            double l = 1;
            double mu = l * u;
            double GaussianNormalDistribution(double t) => 1 / (l * Math.Sqrt(2 * Math.PI)) * Math.Exp(-Math.Pow(t - mu, 2) / 2 * Math.Pow(l, 2));

            LaguerreTabulator Instance = new(GaussianNormalDistribution, 10000, i, j, 15, 0.001);
            (int[] a, double[] b) = Instance.TabulateLaguerreTransformation();
            for (int n = 0; n < a.Length; n++)
            {
                Console.WriteLine($"For lambda={l} u={mu} beta={i} and sigma={j} Transformation value: n={a[n]} transformed=P{b[n]}");
            }

        }
    }
}


public static class Service
{
    public static double Integrate(Func<double, double> f, double a, double b, double margin = 0.001)
    {
        if (a > b)
            throw new ArgumentException("Left bound must be lower than right bound.");

        double area = 0;
        while (Math.Abs(b - a) > margin)
        {
            area += margin * f(a + margin / 2);
            a += margin;
        }
        return area;
    }
}

public class Laguerre
{
    private Func<double, double> _func;
    private double _T;
    private int _beta;
    private int _sigma;
    private int _N;
    private double _epsilon;
    private double[]? _laguerreTransformationValue;
    private double? _experimentValue;

    public Laguerre(Func<double, double> func, double T, int beta, int sigma, int N, double epsilon = 1e-3)
    {
        if (beta > sigma)
            throw new ArgumentException("Sigma must be greater than beta.");

        _func = func;
        this.T = T;
        this.Beta = beta;
        this.Sigma = sigma;
        this.N = N;
        this.Epsilon = epsilon;
    }

    public Func<double, double> Func
    {
        get => _func;
        set => Console.WriteLine("This is a read-only value.");
    }

    public double[] LaguerreTransformationValue
    {
        get
        {
            if (_laguerreTransformationValue is null)
                _laguerreTransformationValue = LaguerreTransformation();
            return _laguerreTransformationValue;
        }
        set => Console.WriteLine("This is a read-only value.");
    }

    public double ExperimentValue
    {
        get
        {
            if (_experimentValue is null)
                _experimentValue = Experiment();
            return _experimentValue.Value;
        }
        set => Console.WriteLine("This is a read-only value.");
    }

    public double T
    {
        get => _T;
        set
        {
            if (value <= 0)
                throw new ArgumentException("T must be positive.");
            Reset();
            _T = value;
        }
    }

    public int Beta
    {
        get => _beta;
        set
        {
            if (value < 0)
                throw new ArgumentException("Beta must be non-negative.");
            Reset();
            _beta = value;
        }
    }

    public int Sigma
    {
        get => _sigma;
        set
        {
            if (value < 0)
                throw new ArgumentException("Sigma must be non-negative.");
            Reset();
            _sigma = value;
        }
    }

    public int N
    {
        get => _N;
        set
        {
            if (value < 1)
                throw new ArgumentException("N must be at least 1.");
            Reset();
            _N = value;
        }
    }

    public double Epsilon
    {
        get => _epsilon;
        set
        {
            if (value <= 0)
                throw new ArgumentException("Epsilon must be positive.");
            Reset();
            _epsilon = value;
        }
    }

    private void Reset()
    {
        _experimentValue = null;
        _laguerreTransformationValue = null;
    }

    public double LaguerreFunction(double t, int n)
    {
        if (_beta < 0 || _beta > _sigma || n < 0)
            throw new ArgumentException("Wrong parameters.");

        double lpp = Math.Sqrt(_sigma) * Math.Exp(-_beta * t / 2);
        double lp = Math.Sqrt(_sigma) * (1 - _sigma * t) * Math.Exp(-_beta * t / 2);

        if (n == 0)
            return lpp;
        if (n == 1)
            return lp;

        for (int i = 2; i <= n; i++)
        {
            double temp = lp;
            lp = ((2 * i - 1 - _sigma * t) * lp / i) - ((i - 1) * lpp / i);
            lpp = temp;
        }

        return lp;
    }

    public double Experiment()
    {
        double[] t = GenerateLinspace(0, T, 1000);
        bool trueForAll;
        double result = 0;

        for (int i = 0; i < t.Length; i++)
        {
            trueForAll = true;
            for (int j = 0; j <= N; j++)
            {
                if (Math.Abs(LaguerreFunction(t[i], j)) >= Epsilon)
                {
                    trueForAll = false;
                    break;
                }
            }

            if (trueForAll)
            {
                result = t[i];
                break;
            }
        }

        if (result == 0)
            throw new Exception("The experiment failed, there is no T that satisfies the condition. You can set a bigger max T to check.");

        _experimentValue = result;
        return result;
    }

    public double[] LaguerreTransformation()
    {
        double to = ExperimentValue;
        double[] result = new double[N + 1];

        for (int k = 0; k <= N; k++)
        {
            result[k] = Service.Integrate(t => _func(t) * LaguerreFunction(t, k) * Math.Exp(-t * (_sigma - _beta)), 0, to, Epsilon);
        }

        _laguerreTransformationValue = result;
        return result;
    }

    public double ReverseLaguerreTransformation(double t)
    {
        double sum = 0;
        for (int k = 0; k < N; k++)
        {
            sum += LaguerreTransformationValue[k] * LaguerreFunction(t, k);
        }
        return sum;
    }

    public static double[] GenerateLinspace(double start, double stop, int num)
    {
        double[] linspace = new double[num];
        double step = (stop - start) / (num - 1);

        for (int i = 0; i < num; i++)
        {
            linspace[i] = start + i * step;
        }

        return linspace;
    }
}

public class LaguerreTabulator : Laguerre
{
    public LaguerreTabulator(Func<double, double> func, double T, int beta, int sigma, int N, double epsilon)
        : base(func, T, beta, sigma, N, epsilon)
    {
        if (beta > sigma)
            throw new ArgumentException("Sigma must be greater than beta.");
    }

    public (double[] t, double[] l) TabulateLaguerre(int n, double to)
    {
        double[] t = GenerateLinspace(1, to, (int)to);
        double[] results = new double[t.Length];

        for (int i = 0; i < t.Length; i++)
        {
            results[i] = LaguerreFunction(t[i], n);
        }

        return (t, results);
    }

    public (double[] t, int[] n, double[] L) TabulateExperiment()
    {
        double t = ExperimentValue;
        double[] ltValues = new double[N + 1];
        for (int n = 0; n <= N; n++)
        {
            ltValues[n] = LaguerreFunction(t, n);
        }

        double[] t_array = new double[N + 1];
        int[] n_array = new int[N + 1];
        for (int i = 0; i <= N; i++)
        {
            t_array[i] = t; n_array[i] = i;
        }

        return (t_array, n_array, ltValues);
    }

    public (int[] n, double[] transformed) TabulateLaguerreTransformation()
    {
        double[] transformed = LaguerreTransformationValue;
        int[] n = new int[N + 1];
        for (int i = 0; i <= N; i++)
        {
            n[i] = i;
        }

        return (n, transformed);
    }
}
