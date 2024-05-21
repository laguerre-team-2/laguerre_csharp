namespace LaguerreClasses;

public static class Service
{
    public static double Integrate(
        Func<double, double> f,
        double a,
        double b,
        double margin = 0.001
    )
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

    public Laguerre(
        Func<double, double> func,
        double T,
        int beta,
        int sigma,
        int N,
        double epsilon = 1e-3
    )
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
            _laguerreTransformationValue ??= LaguerreTransformation();
            return _laguerreTransformationValue;
        }
        set => Console.WriteLine("This is a read-only value.");
    }

    public double ExperimentValue
    {
        get
        {
            _experimentValue ??= Experiment();
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
        double[] t = Service.GenerateLinspace(0, T, 1000);
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
            throw new Exception(
                "The experiment failed, there is no T that satisfies the condition. You can set a bigger max T to check."
            );

        _experimentValue = result;
        return result;
    }

    public double[] LaguerreTransformation()
    {
        double to = ExperimentValue;
        double[] result = new double[N + 1];

        for (int k = 0; k <= N; k++)
        {
            result[k] = Service.Integrate(
                t => _func(t) * LaguerreFunction(t, k) * Math.Exp(-t * (_sigma - _beta)),
                0,
                to,
                Epsilon
            );
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
}

public class LaguerreTabulator : Laguerre
{
    public LaguerreTabulator(
        Func<double, double> func,
        double T,
        int beta,
        int sigma,
        int N,
        double epsilon = 1e-3
    )
        : base(func, T, beta, sigma, N, epsilon)
    {
        if (beta > sigma)
            throw new ArgumentException("Sigma must be greater than beta.");
    }

    public List<(double, double)> TabulateLaguerre(int n, double to, int N = 0)
    {
        List<(double, double)> result = [];
        double[] t;
        if (N <= 0)
        {
            t = Service.GenerateLinspace(0, to, (int)to);
        }
        else
        {
            t = Service.GenerateLinspace(0, to, N);
        }
        Console.WriteLine($"{t[2]}");

        for (int i = 0; i < t.Length; i++)
        {
            result.Add(new(t[i], LaguerreFunction(t[i], n)));
        }

        return result;
    }

    public List<(double, int, double)> TabulateExperiment()
    {
        List<(double, int, double)> result = [];

        for (int i = 0; i <= N; i++)
        {
            result.Add(new(ExperimentValue, i, LaguerreFunction(ExperimentValue, i)));
        }
        return result;
    }

    public List<(int, double)> TabulateLaguerreTransformation()
    {
        List<(int, double)> result = [];
        for (int i = 0; i <= N; i++)
        {
            result.Add(new(i, LaguerreTransformationValue[i]));
        }

        return result;
    }
}
