import java.io.FileNotFoundException;

public class Main
{
    public static void main(String[] args)
    {
        try
        {
            Zadanie.Trapezy();
            Zadanie.RK2();
        }
        catch(FileNotFoundException o)
        {
            System.out.println("File not found");
        }
    }
}