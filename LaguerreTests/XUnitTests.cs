using LaguerreClasses;
using Xunit;

public class AllData : IDisposable
{
    public List<Laguerre> Laguerres { get; }

    public List<LaguerreTabulator> LaguerreTabulators { get; }

    public AllData()
    {
        Laguerres = [];
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

        LaguerreTabulators = [new LaguerreTabulator(func, 10, 2, 5, 10, 0.001)];
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

    [Theory]
    [InlineData(5, 2, 3.97002)]
    public void LaguerreFunction_ValidInput_ReturnsExpectedResult(double t, int n, double expected)
    {
        var laguerre = fixture.Laguerres[0];
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

    [Theory]
    [InlineData(4, -0.001999)]
    [InlineData(1, -0.33491)]
    [InlineData(2, 1.18473)]
    [InlineData(3, 1.29813)]
    [InlineData(5, 0.77927)]
    public void LaguerreTransformation_ValidInput_ReturnsExpectedResult(
        int laguerreIndex,
        double expected
    )
    {
        var laguerre = fixture.Laguerres[laguerreIndex];
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
        Assert.NotNull(result);
    }
}

public class ServiceTests : IClassFixture<AllData>
{
    private readonly AllData fixture;

    public ServiceTests(AllData allData)
    {
        fixture = allData;
    }

    [Theory]
    [InlineData(0, 2, 2.66667)]
    [InlineData(1, 3, 8.66667)]
    [InlineData(-1, 1, 0.66567)]
    public void IntegrationTest_VariousFunctions(double a, double b, double expected)
    {
        Assert.Equal(expected, Service.Integrate((x) => Math.Pow(x, 2), a, b), 5);
    }

    [Theory]
    [InlineData(0, 2, 1.41615)]
    [InlineData(1, 4, 1.19395)]
    public void IntegrationTest_SinFunction(double a, double b, double expected)
    {
        Assert.Equal(expected, Service.Integrate(Math.Sin, a, b), 5);
    }

    [Theory]
    [InlineData(0, 2, 16.0 / 3)]
    [InlineData(-1, 1, 1.33133)]
    public void IntegrationTest_PolynomialFunctions(double a, double b, double expected)
    {
        Assert.Equal(expected, Service.Integrate((x) => 2 * Math.Pow(x, 2), a, b), 5);
    }
}
