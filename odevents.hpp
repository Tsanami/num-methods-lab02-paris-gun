#ifndef ODEVENTS_HPP
#define ODEVENTS_HPP

#include <vector>
#include <functional>
#include <limits>

struct IntegrationEvent {
    enum Type { SURFACE_HIT, MAX_ALTITUDE };
    Type type;
    double time;
    std::vector<double> state;
};

using EventCallback = std::function<void(const IntegrationEvent&)>;

class ODEIntegratorWithEvents {
private:
    EventCallback eventCallback_;
    std::vector<std::vector<std::vector<double>>> trajectoriesBetweenEvents_;
    double surfaceZ_ = 0.0;
    bool checkSurface_ = false;
    bool checkMaxAltitude_ = false;

    double findExactSurfaceTime(const std::vector<double>& before,
                                const std::vector<double>& after) const;

public:
    void addSurfaceEvent(double z_surface = 0.0);
    void addMaxAltitudeEvent();
    void setEventCallback(EventCallback callback);

    std::vector<std::vector<double>> integrateWithEvents(
            const std::vector<double>& initialState,
            double t0, double tEnd,
            double initialStep = 0.01,
            bool useRK4 = true);

    const std::vector<std::vector<std::vector<double>>>& getTrajectoriesBetweenEvents() const;
};

#endif // ODEVENTS_HPP