#include "odevents.hpp"
#include "odesolver.hpp"
#include <algorithm>
#include <iostream>

void ODEIntegratorWithEvents::addSurfaceEvent(double z_surface) {
    checkSurface_ = true;
    surfaceZ_ = z_surface;
}

void ODEIntegratorWithEvents::addMaxAltitudeEvent() {
    checkMaxAltitude_ = true;
}

void ODEIntegratorWithEvents::setEventCallback(EventCallback callback) {
    eventCallback_ = callback;
}

double ODEIntegratorWithEvents::findExactSurfaceTime(const std::vector<double>& before,
                                                     const std::vector<double>& after) const {
    // Линейная интерполяция для нахождения точного времени пересечения поверхности
    const double z_before = before[2] - surfaceZ_;
    const double z_after = after[2] - surfaceZ_;
    return -z_before / (z_after - z_before);
}

std::vector<std::vector<double>> ODEIntegratorWithEvents::integrateWithEvents(
        const std::vector<double>& initialState,
        double t0, double tEnd,
        double initialStep, bool useRK4)
{
    ODESolver solver;
    std::vector<std::vector<double>> currentTrajectory;
    std::vector<double> currentState = initialState;
    double currentTime = t0;
    double currentStep = initialStep;

    double maxAltitude = -std::numeric_limits<double>::max();
    bool maxAltitudeEventFired = false;
    bool surfaceHit = false;

    currentTrajectory.push_back(currentState);

    while (currentTime < tEnd && !surfaceHit) {
        auto newState = useRK4 ?
                        solver.rungeKutta4(currentState, currentTime, currentStep) :
                        solver.euler(currentState, currentTime, currentStep);

        // Проверка на достижение поверхности
        if (checkSurface_ && currentState[2] >= surfaceZ_ && newState[2] < surfaceZ_) {
            double dt = findExactSurfaceTime(currentState, newState);
            auto exactState = useRK4 ?
                              solver.rungeKutta4(currentState, currentTime, dt * currentStep) :
                              solver.euler(currentState, currentTime, dt * currentStep);

            currentTrajectory.push_back(exactState);

            if (eventCallback_) {
                IntegrationEvent event{IntegrationEvent::SURFACE_HIT,
                                       currentTime + dt * currentStep,
                                       exactState};
                eventCallback_(event);
            }

            surfaceHit = true;
            break;
        }

        // Проверка на максимальную высоту
        if (checkMaxAltitude_) {
            if (newState[2] > maxAltitude) {
                maxAltitude = newState[2];
                maxAltitudeEventFired = false;
            } else if (!maxAltitudeEventFired) {
                if (eventCallback_) {
                    IntegrationEvent event{IntegrationEvent::MAX_ALTITUDE,
                                           currentTime,
                                           currentState};
                    eventCallback_(event);
                }
                maxAltitudeEventFired = true;
            }
        }

        currentTrajectory.push_back(newState);
        currentState = newState;
        currentTime += currentStep;

        // Адаптивный выбор шага (упрощенная версия)
        currentStep = std::min(currentStep * 1.1, initialStep);
    }

    trajectoriesBetweenEvents_.push_back(currentTrajectory);
    return currentTrajectory;
}

const std::vector<std::vector<std::vector<double>>>&
ODEIntegratorWithEvents::getTrajectoriesBetweenEvents() const {
    return trajectoriesBetweenEvents_;
}