using LaguerreClasses;

using (StreamWriter writetext = new("data/GaussReverseTransformationData.csv"))
{
    int a = 10;
    int N = 20;
    int beta = 2;
    int sigma = 4;
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

static double func(double x) => 1 / (x + 1);
LaguerreTabulator Instance = new(func, 100, 2, 4, 10, 0.001);
var tabulatedExperimentData = Instance.TabulateExperiment();
using (StreamWriter writetext = new("data/tabulatedExperiment.csv"))
{
    writetext.WriteLine("t,n,L(t)");
    foreach (var value in tabulatedExperimentData)
        writetext.WriteLine($"{value.Item1},{value.Item2},{value.Item3}");
}

var tabulateLaguerreTransformation = Instance.TabulateLaguerreTransformation();
using (StreamWriter writetext = new("data/tabulateLaguerreTransformation.csv"))
{
    writetext.WriteLine("x,t");
    foreach (var value in tabulateLaguerreTransformation)
        writetext.WriteLine($"{value.Item1},{value.Item2}");
}

using (StreamWriter writetext = new("data/lauguerreFuncData.csv"))
{
    int N = 20;
    int a = 3;
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
    List<(string, Laguerre)> laguerres =
    [
        new("funcSonya", new Laguerre(funcSonya, 100, 2, 4, 10)),
        new("funcDemian", new Laguerre(funcDemian, 100, 2, 4, 10)),
        new("funcIvan", new Laguerre(funcIvan, 100, 2, 4, 10)),
        new("funcStefa", new Laguerre(funcStefa, 100, 2, 4, 10)),
        new("funcYuliia", new Laguerre(funcYuliia, 100, 2, 4, 10)),
    ];
    int a = 6;
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
