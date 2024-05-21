using LaguerreClasses;

static double funcSonya(double t)
{
    if (t != 0)
        return Math.Sin(Math.PI / 3) - Math.Atan(t + t / 2) / t;
    else
        return 0;
}

static double funcDemian(double t)
{
    return 2 * (Math.PI / 2 - Math.Atan(t + Math.PI / 6));
}

static double funcIvan(double t)
{
    if (t >= 0 && t <= 2 * Math.PI)
        return Math.PI / 3 - Math.Sin(t + 3 * Math.PI / 2);
    else
        return 0;
}

static double funcStefa(double t)
{
    if (t >= 0 && t <= 2 * Math.PI)
        return -1d / 200d * Math.Sin(t) * Math.Pow(Math.E, t);
    else
        return 0;
}

static double funcYuliia(double t)
{
    return Math.Cos(t) / Math.Pow(t, t);
}

double T = 100;
int beta = 2;
int sigma = 4;
int N = 12;
List<(string, LaguerreTabulator)> LaguerreInstances =
[
    new("funcSonya", new LaguerreTabulator(funcSonya, T, beta, sigma, N)),
    new("funcDemian", new LaguerreTabulator(funcDemian, T, beta, sigma, N)),
    new("funcIvan", new LaguerreTabulator(funcIvan, T, beta, sigma, N)),
    new("funcStefa", new LaguerreTabulator(funcStefa, T, beta, sigma, N)),
    new("funcYuliia", new LaguerreTabulator(funcYuliia, T, beta, sigma, N)),
];

using (StreamWriter writetext = new("data/GaussReverseTransformationData.csv"))
{
    int a = 10;
    List<(string, Laguerre)> laguerres =
    [
        (
            "Gauss(μ=3λ λ=1)",
            new(
                (double t) =>
                    1
                    / (1 * Math.Sqrt(2 * Math.PI))
                    * Math.Exp(-Math.Pow(t - 3, 2) / 2 * Math.Pow(1, 2)),
                100,
                beta,
                sigma,
                N,
                0.001
            )
        ),
        (
            "Gauss(μ=6λ λ=1)",
            new(
                (double t) =>
                    1
                    / (1 * Math.Sqrt(2 * Math.PI))
                    * Math.Exp(-Math.Pow(t - 6, 2) / 2 * Math.Pow(1, 2)),
                100,
                beta,
                sigma,
                N,
                0.001
            )
        ),
        (
            "Gauss(μ=3λ λ=2)",
            new(
                (double t) =>
                    1
                    / (2 * Math.Sqrt(2 * Math.PI))
                    * Math.Exp(-Math.Pow(t - 6, 2) / 2 * Math.Pow(2, 2)),
                100,
                beta,
                sigma,
                N,
                0.001
            )
        ),
        (
            "Gauss(μ=6λ λ=2)",
            new(
                (double t) =>
                    1
                    / (2 * Math.Sqrt(2 * Math.PI))
                    * Math.Exp(-Math.Pow(t - 6, 2) / 2 * Math.Pow(2, 2)),
                100,
                beta,
                sigma,
                N,
                0.001
            )
        ),
    ];
    writetext.Write("x");
    foreach (var item in laguerres)
    {
        writetext.Write($",{item.Item1},{item.Item1}_transformed");
    }
    writetext.WriteLine();
    foreach (double x in Service.GenerateLinspace(1, a, 100 * (a - 1)))
    {
        writetext.Write($"{x}");
        foreach (var laguerre in laguerres)
        {
            writetext.Write(
                $",{laguerre.Item2.Func(x)},{laguerre.Item2.ReverseLaguerreTransformation(x)}"
            );
        }
        writetext.WriteLine();
    }
}

using (StreamWriter writetext = new("data/tabulatedExperiment.csv"))
{
    writetext.Write("n");
    foreach (var item in LaguerreInstances)
    {
        writetext.Write($",t_{item.Item1},L_(t)_{item.Item1}");
    }
    writetext.WriteLine();

    List<List<(double, int, double)>> tabulatedExperiments = [];
    foreach (var laguerre in LaguerreInstances)
    {
        tabulatedExperiments.Add(laguerre.Item2.TabulateExperiment());
    }

    for (int i = 0; i <= N; i++)
    {
        writetext.Write($"{i}");
        foreach (var experiment in tabulatedExperiments)
        {
            writetext.Write($",{experiment[i].Item2},{experiment[i].Item3}");
        }
        writetext.WriteLine();
    }
}

using (StreamWriter writetext = new("data/tabulateLaguerreTransformation.csv"))
{
    writetext.Write("x");
    foreach (var laguerre in LaguerreInstances)
    {
        writetext.Write($",t_{laguerre.Item1}");
    }
    writetext.WriteLine();
    List<List<(int, double)>> tabulateLaguerreTransformations = [];
    foreach (var laguerre in LaguerreInstances)
    {
        tabulateLaguerreTransformations.Add(laguerre.Item2.TabulateLaguerreTransformation());
    }
    for (int i = 0; i <= N; i++)
    {
        writetext.Write($"{i}");
        foreach (var transformation in tabulateLaguerreTransformations)
        {
            writetext.Write($",{transformation[i].Item2}");
        }
        writetext.WriteLine();
    }
}

using (StreamWriter writetext = new("data/lauguerreFuncData.csv"))
{
    int a = 3;
    LaguerreTabulator Instance = new((x) => 1 / (x + 1), 100, 2, 4, 10, 0.001);
    writetext.Write("x");
    for (int i = 0; i < N; i++)
    {
        writetext.Write($",L_{i}");
    }
    writetext.WriteLine();
    foreach (double x in Service.GenerateLinspace(0, a, 1000 * (a - 1)))
    {
        writetext.Write($"{x}");
        for (int n = 0; n < N; n++)
        {
            writetext.Write($",{Instance.LaguerreFunction(x, n)}");
        }
        writetext.WriteLine();
    }
}

using (StreamWriter writetext = new("data/reverseTransformationData.csv"))
{
    int a = 6;
    writetext.Write("x");
    foreach (var item in LaguerreInstances)
    {
        writetext.Write($",{item.Item1},{item.Item1}_transformed");
    }
    writetext.WriteLine();

    foreach (double x in Service.GenerateLinspace(1, a, 100 * (a - 1)))
    {
        writetext.Write($"{x}");
        foreach (var laguerre in LaguerreInstances)
        {
            writetext.Write(
                $",{laguerre.Item2.Func(x)},{laguerre.Item2.ReverseLaguerreTransformation(x)}"
            );
        }
        writetext.WriteLine();
    }
}
