using Xunit;


public class AllData : IDisposable
{
    public List<Laguerre> Laguerres { get; }

    public List<LaguerreTabulator> LaguerreTabulators { get; }

    public AllData()
    {
        Laguerres = new List<Laguerre>();
        static double func(double x) => x * x; // рандомна функ
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
        Laguerres.Add(new Laguerre(func, 10, 2, 5, 10));
        Laguerres.Add(new Laguerre(funcSonya, 100, 2, 4, 10));
        Laguerres.Add(new Laguerre(funcDemian, 100, 2, 4, 10));
        Laguerres.Add(new Laguerre(funcIvan, 100, 2, 4, 10));
        Laguerres.Add(new Laguerre(funcStefa, 100, 2, 4, 10));
        Laguerres.Add(new Laguerre(funcYuliia, 100, 2, 4, 10));

        LaguerreTabulators = new List<LaguerreTabulator>();
        LaguerreTabulators.Add(new LaguerreTabulator(func, 10, 2, 5, 10, 0.001));

    }

    public void Dispose()
    {
        Laguerres.Clear();
        LaguerreTabulators.Clear();
    }
}
public class AllDataTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public AllDataTests(AllData allData)
    {
        fixture = allData;
    }

    [Fact]
    public void Laguerres_NotNull()
    {
        Assert.NotNull(fixture.Laguerres);
    }

    [Fact]
    public void Laguerres_ContainsData()
    {
        Assert.NotEmpty(fixture.Laguerres);
    }

    [Fact]
    public void LaguerreTabulators_NotNull()
    {
        Assert.NotNull(fixture.LaguerreTabulators);
    }

    [Fact]
    public void LaguerreTabulators_ContainsData()
    {
        Assert.NotEmpty(fixture.LaguerreTabulators);
    }
}
public class LaguerreFunctionTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public LaguerreFunctionTests(AllData allData)
    {
        fixture = allData;
    }

    [Fact]
    public void LaguerreFunction_ValidInput_ReturnsExpectedResult()
    {
        var laguerre = fixture.Laguerres[0];
        double t = 5;
        int n = 2;
        double expected = 3.97002;
        double result = laguerre.LaguerreFunction(t, n);
        Assert.Equal(expected, result, 5);
    }

}

public class LaguerreExperimentTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public LaguerreExperimentTests(AllData allData)
    {
        fixture = allData;
    }
    [Fact]
    public void LaguerreExperimentStefaniaFuncTest()
    {
        var laguerre = fixture.Laguerres[4];
        double expected = 43.543543;
        double result = laguerre.ExperimentValue;
        Assert.Equal(expected, result, 5);
    }
    [Fact]
    public void LaguerreExperimentErrorStefaniaFuncTest()
    {
        var laguerre = fixture.Laguerres[4];
        fixture.Laguerres[4].T = 10;
        Action action = () => laguerre.Experiment();
        Assert.Throws<Exception>(action);
    }
}

public class LaguerreTransformationTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public LaguerreTransformationTests(AllData allData)
    {
        fixture = allData;
    }
    [Fact]
    public void LaguerreTransformationStefaniaFuncTest()
    {
        var laguerre = fixture.Laguerres[4];
        double expected = -0.001999;
        double[] result = laguerre.LaguerreTransformationValue;
        Assert.Equal(expected, result[0], 5);
    }
    [Fact]
    public void LaguerreTransformationSonyaFuncTest()
    {
        var laguerre = fixture.Laguerres[1];
        double expected = -0.33491;
        double[] result = laguerre.LaguerreTransformationValue;
        Assert.Equal(expected, result[0], 5);
    }
    [Fact]
    public void LaguerreTransformationDemianFuncTest()
    {
        var laguerre = fixture.Laguerres[2];
        double expected = 1.18473;
        double[] result = laguerre.LaguerreTransformationValue;
        Assert.Equal(expected, result[0], 5);
    }
    [Fact]
    public void LaguerreTransformationIvanFuncTest()
    {
        var laguerre = fixture.Laguerres[3];
        double expected = 1.29813;
        double[] result = laguerre.LaguerreTransformationValue;
        Assert.Equal(expected, result[0], 5);
    }
    [Fact]
    public void LaguerreTransformationYuliaFuncTest()
    {
        var laguerre = fixture.Laguerres[5];
        double expected = 0.77927;
        double[] result = laguerre.LaguerreTransformationValue;
        Assert.Equal(expected, result[0], 5);
    }
}
public class LaguerreTabulatorTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public LaguerreTabulatorTests(AllData allData)
    {
        fixture = allData;
    }

    [Fact]
    public void TabulateLaguerre_ValidInput_ReturnsExpectedResult()
    {
        var tabulator = fixture.LaguerreTabulators.First();
        int n = 2;
        double to = 5;

        var result = tabulator.TabulateLaguerre(n, to);
        Assert.NotNull(result.t);
        Assert.NotNull(result.l);
        Assert.Equal(result.t.Length, result.l.Length);
    }

}
