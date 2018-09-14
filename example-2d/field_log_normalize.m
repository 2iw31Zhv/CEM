function [field_normalized] = field_log_normalize(field, shift)
% 
% 

log_field = log(abs(field)) + shift;
log_field(log_field<0) = 0;

field_normalized = sign(field) .* log_field;
end

