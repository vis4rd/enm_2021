#include <fstream>
#include <concepts>
#include <filesystem>

namespace fu
{
    namespace impl
    {
        template <typename... Pack>
        concept C_ArePrintable = ((std::floating_point<Pack> ||
            std::integral<Pack> ||
            std::convertible_to<Pack, std::string>), ...);
        
        template <typename T>
        concept C_IsArray = requires (const T& t) { t[0]; };

        template <typename T>
        concept C_IsMatrix = requires (const T& t) { t[0][0]; };

        template <typename... Pack>
        concept C_AreArray = (C_IsArray<Pack> && ...);

        template <typename... Pack>
        concept C_AreMatrix = (C_IsMatrix<Pack> && ...);

        template <typename T>
        requires C_IsArray<T>
        constexpr auto ArrayLength(const T& arr)
        {
            return (sizeof(arr)/sizeof(arr[0]));
        }

        template <typename First, typename... Rest>
        requires C_IsArray<First> && C_AreArray<Rest...>
        constexpr auto ArrayLengthInPack(First&& arr, [[maybe_unused]] Rest&& ...rest)
        {
            return ArrayLength(std::forward<First>(arr));
        }

        template <typename T>
        using RR = std::remove_reference_t<T>;

        template <typename Stream, typename... Pack>
        requires std::derived_from<RR<Stream>, std::basic_ios<typename Stream::char_type>>
        constexpr bool writeToStream(Stream& strout, Pack&& ...vars)
        {
            ((strout << vars << " "), ...);
            strout << std::endl;

            return true;
        }

        template <typename Fstream, typename... Pack>
        requires std::same_as<RR<Fstream>, std::ofstream>
        constexpr bool appendToFile(Fstream&& fstrout, Pack&& ...vars)
        {
            impl::writeToStream(fstrout, std::forward<Pack>(vars)...);
            return true;
        }
    }  // namespace impl

    template <typename... Pack>
    requires impl::C_AreArray<Pack...>
    constexpr bool saveArraysToFile(const std::string& filename, Pack&& ...vars)
    {
        std::ofstream file(filename);
        for(long unsigned int i{0}; i < impl::ArrayLengthInPack(std::forward<Pack>(vars)...); i++)
        {
            impl::appendToFile(file, std::forward<Pack>(vars)[i]...);
        }
        file.close();
        return true;
    }

    template <typename... Pack>
    requires impl::C_AreMatrix<Pack...>
    constexpr bool saveMatricesToFile(const std::string& filename, Pack&& ...vars)
    {
        std::ofstream file(filename);
        for(long unsigned int i{0}; i < impl::ArrayLengthInPack(std::forward<Pack>(vars)...); i++)
        {
            for(long unsigned int j{0}; j < impl::ArrayLengthInPack(std::forward<Pack>(vars)...); j++)
            {
                impl::appendToFile(file, std::forward<Pack>(vars)[i][j]...);
            }
        }
        file.close();
        return true;
    }

    template <typename... Pack>
    requires impl::C_AreArray<Pack...>
    constexpr bool AreSameLength(Pack&& ...vars)
    {
        const auto sum = (impl::ArrayLength(std::forward<Pack>(vars)) + ...);
        const auto count = sizeof...(vars);
        const auto res = (sum/count);
        return ((res == impl::ArrayLength(std::forward<Pack>(vars))) && ...);
    }

    template <typename T>
    requires std::same_as<impl::RR<T>, std::string> ||
        std::convertible_to<impl::RR<T>, std::string>
    constexpr bool removeFile(T&& filename)
    {
        return std::filesystem::remove(filename);
    }
}
