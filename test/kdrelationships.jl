linear_response = response(Linear())

@test linear_response(0.0, 1.0) == 0.0
@test linear_response(1.0, 1.0) == 1.0
@test linear_response(0.0, -1.0) == 0.0
@test linear_response(1.0, -1.0) == -1.0

sigmoidal_response = response(Sigmoidal())

@test sigmoidal_response(0.0, 1.0) == 0.0
@test sigmoidal_response(1.0, 1.0) == 1.0
@test sigmoidal_response(0.0, -1.0) == -0.0
@test sigmoidal_response(1.0, -1.0) == -1.0
